import jawm
import os
from pathlib import Path

AGEPY_IMAGE="mpgagebioinformatics/agepy:9acb5de"

annotations=jawm.Process(
    name="annotations",
    when=lambda p: not os.path.isfile( os.path.join( p.var["deseq2_output"],"annotated",  "biotypes_go.txt") ) ,
    script="""\
#!/usr/local/bin/python
import pandas as pd
import os
import AGEpy as age
import shutil

if not os.path.isdir("{{deseq2_output}}/annotated/") :
    os.makedirs("{{deseq2_output}}/annotated/")

if not os.path.isfile("{{biotypes_go}}"):
    print("Not found: {{biotypes_go}}")
    print("Contacting the BiomartServer.")
    from biomart import BiomartServer
    attributes=["ensembl_gene_id","external_gene_name","go_id","name_1006"]
    
    if "{{biomart_host}}" :
        biomart_host="{{biomart_host}}"
    else:
        biomart_host=age.get_ensembl_biomart_archive_url("{{release}}")
        
    print(f"biomart host: {biomart_host}" )
    server = BiomartServer( biomart_host )

    if "{{biomart_dataset}}" :
        biomart_dataset="{{biomart_dataset}}"
    elif "{{organism}}" :
        organism="{{organism}}"
        biomart_dataset=organism.split("_")[0][0] + organism.split("_")[1] + "_gene_ensembl"
    else:
        raise Exception("You must provide a biomart_dataset or an organism.")

    organism=server.datasets[ biomart_dataset ]
    response=organism.search({"attributes":attributes})
    response=response.content.decode().split("\\n")
    response=[s.split("\\t") for s in response ]
    bio_go=pd.DataFrame(response,columns=attributes)
    bio_go=bio_go.sort_values(by=attributes, ascending=True )
    bio_go.to_csv("{{biotypes_go}}".replace("biotypes_go.txt", "biotypes_go_raw_topgo.txt"), index=None, sep="\\t")
    bio_go.to_csv("{{deseq2_output}}/annotated/biotypes_go_raw_topgo.txt", index=None, sep="\\t")
    bio_go=bio_go[["ensembl_gene_id","go_id","name_1006"]]
    bio_go.to_csv("{{biotypes_go}}".replace("biotypes_go.txt", "biotypes_go_raw.txt"), index=None, sep="\\t")
    bio_go.to_csv("{{deseq2_output}}/annotated/biotypes_go_raw.txt", index=None, sep="\\t")
    bio_go.columns = ["ensembl_gene_id","GO_id","GO_term"]
    bio_go=bio_go.sort_values(by=["ensembl_gene_id","GO_id"], ascending=True )

    def join_sorted_unique(series):
        # drop NaNs, cast to str once, dedup while preserving first-seen order, then sort for full determinism
        vals = pd.Series(series).dropna().astype(str).unique()   # preserves first appearance
        # vals = sorted(vals)                                      # or omit this line if you prefer original order
        return "; ".join(vals)

    def CombineAnn(df):
        return pd.Series(dict(ensembl_gene_id = ("ensembl_gene_id", join_sorted_unique) ,\
                        GO_id = ("GO_id", join_sorted_unique) ,\
                        GO_term = ("GO_term", join_sorted_unique)  ,\
                        ) )

    bio_go = (
        bio_go.groupby("ensembl_gene_id", sort=False)
            .agg(
                ensembl_gene_id=("ensembl_gene_id", "first"),
                GO_id=("GO_id", join_sorted_unique),
                GO_term=("GO_term", join_sorted_unique),
            )
            .reset_index(drop=True)
            # ensure column order
            [["ensembl_gene_id","GO_id","GO_term"]]
    )

    bio_go=bio_go.dropna(subset=["ensembl_gene_id"])

    GTF=age.readGTF("{{gtf}}")
    GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
    GTF["gene_biotype"]=age.retrieve_GTF_field(field="gene_biotype",gtf=GTF)
    GTF=GTF[["gene_id","gene_biotype"]].drop_duplicates().dropna(subset=["gene_id"])
    GTF.columns=["ensembl_gene_id","gene_biotype"]
    GTF=GTF.astype(str)
    bio_go=bio_go.astype(str)
    bio_go=pd.merge(GTF,bio_go,on=["ensembl_gene_id"],how="outer")

    # except Exception as e:
    #     print("An error occurred:", e,"\\ncontinuing without biomart go annotations.")
    #     bio_go=pd.DataFrame(columns=["ensembl_gene_id"])

    bio_go=bio_go.sort_values(by=["ensembl_gene_id"],ascending=True)
    bio_go.to_csv("{{deseq2_output}}/annotated/biotypes_go.txt", sep= "\\t", index=None)
    bio_go.to_csv("{{biotypes_go}}", sep= "\\t", index=None)     

else:
    shutil.copy( "{{biotypes_go}}".replace("biotypes_go.txt", "biotypes_go_raw_topgo.txt"), "{{deseq2_output}}/annotated/biotypes_go_raw_topgo.txt" )
    shutil.copy( "{{biotypes_go}}".replace("biotypes_go.txt", "biotypes_go_raw.txt"), "{{deseq2_output}}/annotated/biotypes_go_raw.txt" )
    shutil.copy( "{{biotypes_go}}", "{{deseq2_output}}/annotated/biotypes_go.txt")
""",
    var={
        "biomart_dataset": "",
        "organism":"",
        "biomart_host":""
    },
    desc={
        "gtf":"",
        "biomart_dataset": "",
        "release": "",
        "deseq2_output":"",
        "biotypes_go":"'<path/to/biotypes_go.txt' inclusive 'biotypes_go.txt'"
    },
    container=AGEPY_IMAGE,
    manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "1:00:00" }
)


parse_submission=jawm.Process(
    name="parse_submission",
    when=lambda p: not os.path.isfile( os.path.join( p.var["deseq2_output"], "models.txt") ) ,
    script="""\
#!/usr/local/bin/python
import pandas as pd
from biomart import BiomartServer
import itertools
import AGEpy as age
import os

sfile="{{samplestable}}"
if os.path.isfile(sfile):
    sdf=pd.read_excel(sfile)

elif os.path.isdir("{{kallisto_output}}"):
    print("Could not find a samples table of the form `pd.DataFrame(columns=['Files/Folders','group'])`")
    print("Working on predefined / expected kallisto output to generate test comparisons.")
    folders=os.listdir("{{kallisto_output}}")
    folders=[ s for s in folders if os.path.isdir(f"{{kallisto_output}}/{s}") ]
    folders=[ s for s in folders if not s.startswith("tmp.") ]
    folders=[ s for s in folders if "{{rep_prefix}}" in s ]
    if not folders :
        raise Exception("Could not find a samples table of the form `pd.DataFrame(columns=['Files/Folders','group'])` nor a kallisto folder with Replicate information eg. <group>.Rep_<number> in it's sub folder names")
    folders.sort()
    groups=[ s.split("{{rep_prefix}}")[0] for s in folders ]
    sdf=pd.DataFrame( { "Files/Folders":folders , "group":groups } )

sam_df=sdf.copy()
files_col=sam_df.columns.tolist()[0]
cond_col=sam_df.columns.tolist()[1]
sam_df.index=[s.split("{{read1_suffix}}")[0] for s in sam_df[files_col].tolist()]
sam_df=sam_df[[cond_col]]
sam_df.to_csv("{{deseq2_output}}/samples_MasterTable.txt", sep="\\t")
fs=sdf.columns.tolist()[0]
sdf[fs]=sdf[fs].apply(lambda x: x.split("{{read1_suffix}}")[0] )
sdf.index=sdf[fs].tolist()
sdf=sdf.drop([fs],axis=1)
cols=sdf.columns.tolist()
mods=[ [x,y] for x in cols for y in cols ]
interactions=[]
for c in mods:
    if c not in interactions:
        if [c[1], c[0]] not in interactions:
            if c[0] != c[1]:
                interactions.append(c)
single_models=cols
textout=[]
for m in single_models:
    variants=[ s for s in cols if s != m ] 
    tmp=sdf.copy()
    tmp["_group_"]=m
    for c in variants:
        tmp["_group_"]=tmp["_group_"].astype(str)+"."+c+"_"+tmp[c].astype(str)
    for g in list(set(tmp["_group_"].tolist())):
        outdf=tmp[tmp["_group_"]==g]
        print(outdf)
        model_data=list(set(outdf[m].tolist()))
        model_pairs=[ list(set([x,y])) for x in model_data for y in model_data ]
        model_pairs_=[ pair for pair in model_pairs if len(pair) > 1 ]
        model_pairs=[]
        for pair in model_pairs_:
            if pair not in model_pairs:
                model_pairs.append(pair)
        for pair in model_pairs:
            outdf_=outdf[outdf[m].isin(pair)]
            outdf_=outdf_.drop(["_group_"],axis=1)
            coef=outdf_[m].tolist()
            ref=coef[0]
            target=[ t for t in coef if t != ref ][0]
            coef=m+"_"+str(target)+"_vs_"+str(ref)
            filename=coef+g.split(m)[-1]
            print("{{deseq2_output}}/"+filename+".input.tsv")
            outdf_.to_csv("{{deseq2_output}}/"+filename+".input.tsv", sep="\t")
            text=[ filename, m, g, str(ref), coef ]
            text="\\t".join(text)
            textout.append(text)
        
with open("{{deseq2_output}}/models.txt", "w") as mout:
    mout.write("\\n".join(textout) + "\\n")
""",
    var={
        "samplestable": "",
        "rep_prefix":".Rep_"
    },
    desc={
        "samplestable":"",
        "deseq2_output":"",
        "read1_suffix":""
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "4GB", "-t": "1:00:00" }
)

tx2gene=jawm.Process(
    name="tx2gene",
    when=lambda p: not os.path.isfile( os.path.join( p.var["deseq2_output"] , "tx2gene.csv") ) ,
    script="""\
#!/usr/local/bin/python
import AGEpy as age
import os
GTF=age.readGTF("{{gtf}}")
GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
GTF["transcript_id"]=age.retrieve_GTF_field(field="transcript_id",gtf=GTF)
tx2gene=GTF[["transcript_id","gene_id"]].drop_duplicates().dropna()
tx2gene.columns=["TXNAME","GENEID"]
tx2gene[["TXNAME","GENEID"]].to_csv("{{deseq2_output}}/tx2gene.csv", quoting=1, index=None)
""",
    desc={
        "gtf":"",
        "deseq2_output":"", 
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "4GB", "-t": "1:00:00" }
)


tx2gene_proc=jawm.Process(
    name="tx2gene_proc",
    when=lambda p: not os.path.isfile( os.path.join( p.var["deseq2_output"] , "deseq2.part1.Rdata") ) ,
    script="""\
#!/usr/bin/Rscript
library(tidyverse)
library(tximportData)
library(tximport)
library(DESeq2)
library(rhdf5)
library(readr)
library(apeglm)
tx2gene <- read_csv("{{deseq2_output}}/tx2gene.csv")
tx2gene <- tx2gene[rowSums(is.na(tx2gene)) == 0,]
tx2gene <- droplevels(tx2gene)
save.image("{{deseq2_output}}/deseq2.part1.Rdata")
sessionInfo()
""",
    desc={
        "project_folder":""
    },
    container="mpgagebioinformatics/deseq2:1.38.0",
    manager_slurm={ "-c": 1, "--mem": "10GB", "-t": "1:00:00" }
)

deseq2=jawm.Process(
    name="deseq2",
    when=lambda p: not os.path.isfile( os.path.join( p.var["deseq2_output"] , p.var["input_file"].split(".input.tsv")[0]+".results.tsv" ) ) ,
    script="""\
#!/usr/bin/Rscript
load("{{deseq2_output}}/deseq2.part1.Rdata")
library(tximportData)
library(tximport)
library(DESeq2)
library(rhdf5)
library(readr)
library(apeglm)

print("{{input_file}}.tsv")

ref=stringr::str_split("{{input_file}}", "_vs_")[[1]][[2]]
ref=stringr::str_split(ref, ".input.tsv")[[1]][[1]]

coef=stringr::str_split("{{input_file}}", ".input.tsv")[[1]][[1]]

out=paste("{{deseq2_output}}/",coef,".results.tsv", sep="")
out_vsd=paste("{{deseq2_output}}/",coef,".vsd.counts.tsv", sep="")

# filein
sampleTable<-read.delim2("{{deseq2_output}}/{{input_file}}",sep = "\t", row.names = 1)
samples<-row.names(sampleTable)

# dir
dir<-"{{kallisto_output}}"
files<-file.path(dir, samples, "abundance.tsv")
print(files)
names(files)<-samples
txi <- tximport(files, type = "kallisto", txOut = TRUE)
txi <- summarizeToGene(txi, tx2gene = tx2gene)

# model
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~group )

# if circRNA file is present, add circRNA counts to dds, else proceed
# here, get count table,
# add circRNA counts
circRNA_folder="{{circRNA}}"
if(circRNA_folder != "None"){
  # gene count table
  gene_count = as.data.frame(counts(dds, normalized=FALSE))
  # circRNA table
  circRNA = read.delim('{{deseq2_output}}/{{circRNA}}/CircRNACount', check.names = FALSE, as.is = TRUE )
  coord = read.delim('{{deseq2_output}}/{{circRNA}}/CircCoordinates', check.names = FALSE, as.is = TRUE)
  # reformat circRNA header
  names(circRNA) <- gsub('.Chimeric.out.junction', '', names(circRNA))
  row.names(circRNA) <- paste0('circ_', coord[, 'Chr'], ':', coord[,'Start'], '-', coord[,'End'], '|', coord[,'Strand'], '|', coord[,'Gene'])
  circRNA = circRNA[,names(gene_count)]
  # add circRNAs to gene count table
  gene_count = rbind(gene_count, circRNA)
  # create DESeqDataSetFromMatrix
  dds <- DESeqDataSetFromMatrix(gene_count, sampleTable, ~group)
}

# model
dds$group <- relevel(dds$group, ref = ref)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res.counts<-counts(dds, normalized=TRUE)


## save vst counts
vsd <- vst(dds, blind = TRUE)
vsd <- assay(vsd)
vsd <- as.data.frame(vsd)
vsd$ensembl_gene_id <- rownames(vsd)
vsd <- vsd[, c("ensembl_gene_id", setdiff(colnames(vsd), "ensembl_gene_id"))]
write.table( vsd, file = out_vsd, sep = "\\t", quote = FALSE, row.names = FALSE )


# coef
resLFC <- lfcShrink(dds,  coef=coef, type="apeglm")
counts.dge<-merge(res.counts,resLFC, by=0, all=FALSE)
row.names(counts.dge)<-counts.dge$Row.names
counts.dge<-counts.dge[, !(colnames(counts.dge) %in% c("Row.names"))]
counts.dge <- counts.dge[order(counts.dge$pvalue),]

# file out
write.table(counts.dge, out, sep="\\t")

sessionInfo()
""",
    var={
        "circRNA":"None"
    },
    desc={
        "input_file":"",
        "deseq2_output":"",
        "kallisto_output":"",
        "circRNA": ""
    },
    container="mpgagebioinformatics/deseq2:1.38.0",
    manager_slurm={ "-c": 20, "--mem": "20GB", "-t": "1:00:00" }
)

mastertable=jawm.Process(
    name="mastertable",
    when=lambda p: not os.path.isfile( os.path.join( p.var["deseq2_output"] , "all_results_stats.xlsx" ) ) ,
    script="""\
#!/usr/bin/Rscript
load("{{deseq2_output}}/deseq2.part1.Rdata")
library(tidyverse)
library(tximportData)
library(tximport)
library(DESeq2)
library(rhdf5)
library(readr)
library(apeglm)

sampleTable<-read.delim2("{{deseq2_output}}/samples_MasterTable.txt",sep = "\\t", row.names = 1)
samples<-row.names(sampleTable)
# dir
dir <- "{{kallisto_output}}"
files<-file.path(dir, samples, "abundance.tsv")
names(files)<-samples
txi <- tximport(files, type = "kallisto", txOut = TRUE)
txi <- summarizeToGene(txi, tx2gene = tx2gene)
tpm <- txi$abundance

# model
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ group )
# if circRNA file is present, add circRNA counts to dds, else proceed
# here, get count table, 
# add circRNA counts
circRNA_folder="{{circRNA}}"
if(circRNA_folder != "None"){
  # gene count table
  gene_count = as.data.frame(counts(dds, normalized=FALSE))
  # circRNA table
  circRNA = read.delim('{{deseq2_output}}/{{circRNA}}/CircRNACount', check.names = FALSE, as.is = TRUE )
  coord = read.delim('{{deseq2_output}}/{{circRNA}}/CircCoordinates', check.names = FALSE, as.is = TRUE)

  # reformat circRNA header
  names(circRNA) <- gsub('.Chimeric.out.junction', '', names(circRNA))
  row.names(circRNA) <- paste0('circ_', coord[, 'Chr'], ':', coord[,'Start'], '-', coord[,'End'], '|', coord[,'Strand'], '|', coord[,'Gene'])

  circRNA = circRNA[,names(gene_count)]

  # add circRNAs to gene count table
  gene_count = rbind(gene_count, circRNA)
  # create DESeqDataSetFromMatrix
  dds <- DESeqDataSetFromMatrix(gene_count, sampleTable, ~ group)
} 
dds <- estimateSizeFactors(dds)
res.counts <- counts(dds, normalized=TRUE)
res_counts <- as.data.frame(res.counts)
openxlsx::write.xlsx(res_counts, "{{deseq2_output}}/all_res_counts.xlsx", row.names = TRUE, col.names = TRUE)
write.table(res_counts, "{{deseq2_output}}/all_res_counts.tsv", sep = "\\t", quote = T, row.names = T)
result_tables <- list.files('{{deseq2_output}}/', pattern = '.results.tsv')
for(f in result_tables){
  tmp <- read.delim(paste0('{{deseq2_output}}/', f))
  tmp <- tmp[,c('log2FoldChange', 'pvalue', 'padj')]
  names(tmp) <- paste(names(tmp), gsub('.results.tsv', '', gsub('Group_', '', f)), sep = '.')
  res_counts <- merge(res_counts, tmp, by.x = 'row.names', by.y = 'row.names', all = TRUE)
  #print(nrow(res_counts))
  names(res_counts)[names(res_counts) == "Row.names"] <- "ensembl_gene_id"
  print(length(res_counts[,'ensembl_gene_id']))
  row.names(res_counts) <- res_counts[,'ensembl_gene_id']
  res_counts <- res_counts[, !duplicated(colnames(res_counts))]
}

# calculate fpkm values
fpkm.deseq = as.data.frame(fpkm(dds))
names(fpkm.deseq) = paste0('fpkm.', names(fpkm.deseq))
res_counts = merge(res_counts, fpkm.deseq, by = 'row.names', all = TRUE)
# Ensure ensembl_gene_id column exists and is first column

if (!"ensembl_gene_id" %in% colnames(res_counts)) {
    if ("Row.names" %in% colnames(res_counts)) {
        res_counts$ensembl_gene_id <- res_counts$Row.names
    } else {
        # fallback: use actual rownames if Row.names column is not present
        res_counts$ensembl_gene_id <- rownames(res_counts)
    }
}

if ("Row.names" %in% colnames(res_counts)) {
  res_counts <- res_counts[, !(colnames(res_counts) %in% "Row.names")]
}
res_counts <- res_counts[, c("ensembl_gene_id", setdiff(colnames(res_counts), "ensembl_gene_id"))]

write.table(res_counts, "{{deseq2_output}}/all_results_stats.tsv", sep = "\\t", quote = F, row.names = F)
openxlsx::write.xlsx(res_counts, "{{deseq2_output}}/all_results_stats.xlsx", row.names = F, col.names = TRUE)


# Convert to data frame and keep gene IDs
tpm_df <- as.data.frame(tpm)
tpm_df$ensembl_gene_id <- rownames(tpm_df)

# Move gene_id to first column
tpm_df <- tpm_df[, c("ensembl_gene_id", setdiff(colnames(tpm_df), "ensembl_gene_id"))]

write.table(tpm_df, "{{deseq2_output}}/tpm.tsv", sep = "\\t", quote = F, row.names = F)


## save vst counts
vsd <- vst(dds, blind = TRUE)
vsd <- assay(vsd)
vsd <- as.data.frame(vsd)
vsd$ensembl_gene_id <- rownames(vsd)
vsd <- vsd[, c("ensembl_gene_id", setdiff(colnames(vsd), "ensembl_gene_id"))]
write.table( vsd, file = "{{deseq2_output}}/all.samples.vsd.counts.tsv", sep = "\\t", quote = FALSE, row.names = FALSE )

## raw counts
raw_counts <- counts(dds, normalized = FALSE)
df_raw <- as.data.frame(raw_counts)
# Add gene IDs
df_raw$ensembl_gene_id <- rownames(df_raw)
# Force it to be the first column
df_raw <- df_raw[, c("ensembl_gene_id", setdiff(colnames(df_raw), "ensembl_gene_id"))]
# Write to TSV
write.table( df_raw, file = "{{deseq2_output}}/all.samples.raw.counts.tsv", sep = "\\t", quote = FALSE, row.names = FALSE )

sessionInfo()
""",
    var={
        "circRNA":"None"
    },
    desc={
        "circRNA": "",
        "deseq2_output":"",
        "kallisto_output":"",
    },
    container="mpgagebioinformatics/deseq2:1.38.0",
    manager_slurm={ "-c": 8, "--mem": "40GB", "-t": "4:00:00" }
)

annotator=jawm.Process(
    name="annotator",
    when=lambda p: not os.path.isfile( os.path.join( p.var["deseq2_output"], "annotated", "masterTable_annotated.xlsx" ) ) ,
    script="""\
#!/usr/local/bin/python
import pandas as pd
import os
import AGEpy as age
import shutil

if not os.path.isdir("{{deseq2_output}}/annotated/") :
    os.makedirs("{{deseq2_output}}/annotated/")

if os.path.exists("{{deseq2_output}}/annotated/biotypes_go.txt"):
    bio_go=pd.read_csv( "{{deseq2_output}}/annotated/biotypes_go.txt", sep="\t")
else:
    bio_go=pd.DataFrame(columns=["ensembl_gene_id"])

GTF=age.readGTF("{{gtf}}")
GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
GTF["gene_name"]=age.retrieve_GTF_field(field="gene_name",gtf=GTF)
id_name=GTF[["gene_id","gene_name"]].drop_duplicates()
id_name.reset_index(inplace=True, drop=True)
id_name.columns=["ensembl_gene_id","gene_name"]

vsd_files=os.listdir("{{deseq2_output}}/")
vsd_files=[ s for s in vsd_files if "vsd.counts.tsv" in s ]
for f in vsd_files:
    df=pd.read_csv("{{deseq2_output}}/"+f,sep="\\t")
    df=pd.merge( id_name, df, on=["ensembl_gene_id"], how="right" ) # change to gene_id
    df.to_excel( "{{deseq2_output}}/annotated/" + f.replace('.tsv', '.xlsx') , index=None )

deg_files=os.listdir("{{deseq2_output}}/")
deg_files=[ s for s in deg_files if "results.tsv" in s ]
i=1
s=[]
dfs={}
for f in deg_files:
    df=pd.read_table("{{deseq2_output}}/"+f)
    df=pd.merge(id_name,df,left_on=["ensembl_gene_id"],right_index=True, how="right") # change to gene_id
    df=pd.merge(df,bio_go,on=["ensembl_gene_id"],how="left")
    df=df.sort_values(by=["padj"],ascending=True)
    df.to_csv("{{deseq2_output}}/annotated/"+f, sep="\\t",index=None)
    df.to_excel("{{deseq2_output}}/annotated/"+f.replace('.tsv', '.xlsx'), index=None)
    n=f.split(".results.tsv")[0]
    s.append([i,n])
    df=df[df["padj"]<0.05]
    df.reset_index(inplace=True, drop=True)
    dfs[i]=df
    i=i+1
sdf=pd.DataFrame(s,columns=["sheet","comparison"])
EXC=pd.ExcelWriter("{{deseq2_output}}/annotated/significant.xlsx")
sdf.to_excel(EXC,"summary",index=None)
for k in list(dfs.keys()):
    dfs[k].to_excel(EXC, str(k),index=None)
EXC.close()
mt=pd.read_csv("{{deseq2_output}}/all_results_stats.tsv", sep="\\t")
mt_ann=pd.merge(id_name,mt,on=["ensembl_gene_id"], how="right")
mt_ann=pd.merge(mt_ann,bio_go,on=["ensembl_gene_id"],how="left")
mt_ann.to_csv("{{deseq2_output}}/annotated/masterTable_annotated.tsv", sep="\\t",index=None)
mt_ann.to_excel("{{deseq2_output}}/annotated/masterTable_annotated.xlsx", index=None)

tpm=pd.read_csv("{{deseq2_output}}/tpm.tsv", sep="\\t")
tpm=pd.merge(id_name,tpm,on=["ensembl_gene_id"], how="right")
tpm.to_excel("{{deseq2_output}}/annotated/tpm.xlsx", index=None)

""",
    desc={
        "gtf":"",
        "biomart_host": "",
        "deseq2_output":"",
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 2, "--mem": "4GB", "-t": "1:00:00" }
)

david=jawm.Process(
    name="david",
    when=lambda p: (  ( not os.path.isfile( os.path.join( p.var["deseq2_output"] , "annotated", "david.touch" ) ) )  &  ( p.var["DAVIDUSER"] != "" ) ) ,
    script="""\
#!/usr/local/bin/python
import pandas as pd
import AGEpy as age
import os 
import sys
from pathlib import Path

file_path = Path("{{deseq2_output}}/annotated/david.touch")

deseq2="{{deseq2_output}}/annotated/"
files=os.listdir(deseq2)
files=[ s for s in files if ".results.tsv" in s ]
for f in files:
    if os.path.isfile(deseq2+f.replace("results.tsv","DAVID.xlsx")):
        continue
    df=pd.read_csv(deseq2+f,sep="\\t")
    df=df[df["padj"]<0.05]
    if df.shape[0] > 0:
        dics=df[["ensembl_gene_id","gene_name","log2FoldChange"]]
        dics["ensembl_gene_id"]=dics["ensembl_gene_id"].apply(lambda x: x.upper())
        dics.index=dics["ensembl_gene_id"].tolist()
        names_dic=dics[["gene_name"]].to_dict()["gene_name"]
        exp_dic=dics[["log2FoldChange"]].to_dict()["log2FoldChange"]

        genes=df["ensembl_gene_id"].tolist()

        try:

            DAVID=age.DAVIDenrich(database="{{daviddatabase}}",\
                        categories="GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,KEGG_PATHWAY,PFAM,PROSITE,GENETIC_ASSOCIATION_DB_DISEASE,OMIM_DISEASE",\
                    user="{{DAVIDUSER}}",\
                    ids=genes, verbose=True)

            if os.path.exists( deseq2+f.replace("results.tsv","DAVID.failed") ) :
                os.remove( deseq2+f.replace("results.tsv","DAVID.failed") )

        except:
            DAVID=""
            print( f.replace("results.tsv","DAVID.failed") )
            david_file_path = Path( deseq2+f.replace("results.tsv","DAVID.failed") )

        if type(DAVID) == type(pd.DataFrame()):
            #for c in DAVID.columns.tolist():
            #    DAVID[c]=DAVID[c].apply(lambda x: x.decode())
            #print(DAVID.head(),DAVID["geneIds"].tolist(),names_dic )
            DAVID["genes name"]=DAVID["geneIds"].apply(lambda x: ", ".join([ str(names_dic[s.upper()]) for s in x.split(", ") ] ) )
            DAVID["log2fc"]=DAVID["geneIds"].apply(lambda x: ", ".join([ str(exp_dic[s.upper()]) for s in x.split(", ") ] ) )

            DAVID.to_csv(deseq2+f.replace("results.tsv","DAVID.tsv"), sep="\\t", index=None)
            EXC=pd.ExcelWriter(deseq2+f.replace("results.tsv","DAVID.xlsx"))
            df.to_excel(EXC,"genes",index=None)
            for cat in DAVID["categoryName"].tolist():
                tmp=DAVID[DAVID["categoryName"]==cat]
                tmp.to_excel(EXC,cat,index=None)
            EXC.close()

file_path.touch()
""",
    var={
        "DAVIDUSER":""
    },
    desc={
        "DAVIDUSER": "",
        "daviddatabase": "",
        "deseq2_output": ""
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "12:00:00" }
)

topgo=jawm.Process(
    name="topgo",
    when=lambda p: not os.path.isfile( os.path.join( p.var["deseq2_output"] , "annotated", p.var["input_file"].replace("results.tsv","topGO.tsv") ) ) ,
    script="""\
#!/usr/bin/Rscript
library(topGO)
library(biomaRt)
library(plyr)
library(openxlsx)
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fin="{{deseq2_output}}/annotated/{{input_file}}" # topGO.tsv
Din = read.delim(fin, as.is = TRUE)
# head(Din)
# make geneList input for topGO
sigGenes = subset(Din, padj <= 0.05)[, 'ensembl_gene_id']
geneList = factor(as.integer(Din[,'ensembl_gene_id'] %in% sigGenes))

# print( "pipeline only go through when there is significant genes")
if ( length(levels(geneList)) > 1) {

    names(geneList) = Din[,'ensembl_gene_id']
    results = list(genes = subset(Din, padj <= 0.05))
    ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # This part if you are working in R version 4 or higher
    # biomartCacheClear()
    # ensembl <- useEnsembl(biomart = "genes", dataset = "{{biomart_dataset}}", host="{{biomart_host}}")
    # host=stringr::str_split("{{biomart_host}}", "/biomart")[[1]][[1]]
    # ensembl = useMart("ensembl",dataset="{{biomart_dataset}}", host=host)
    # goterms = getBM(attributes = c("ensembl_gene_id", 
    #                                        "external_gene_name",
    #                                        "go_id", "name_1006"),  
    #                          filters = 'ensembl_gene_id',
    #                          values = Din[, 'ensembl_gene_id'],
    #                          mart = ensembl)

    # print("goterms<-read.delim")
    goterms<-read.delim("{{deseq2_output}}/annotated/biotypes_go_raw_topgo.txt", as.is = TRUE)
    ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    geneID2GO = list()
    for(i in 1:nrow(goterms)){
        if(goterms[i, 'ensembl_gene_id'] %in% names(geneID2GO)){
            geneID2GO[[goterms[i, 'ensembl_gene_id']]] = c(geneID2GO[[goterms[i, 'ensembl_gene_id']]], 
                                                                    goterms[i, 'go_id'])
        } else {
            geneID2GO[[goterms[i, 'ensembl_gene_id']]] = c(goterms[i, 'go_id'])
        }
    }

    GOTERM_categories = list('BP' = 'GOTERM_BP',
                            'MF' = 'GOTERM_MF', 
                            'CC' = 'GOTERM_CC')

    for(go_cat in names(GOTERM_categories)){
        ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        sampleGOdata <- new("topGOdata",
                            description = "Simple session",
                            ontology = go_cat,
                            allGenes = geneList,
                            nodeSize = 10,
                            annot = annFUN.gene2GO, 
                            gene2GO = geneID2GO)


        ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")


        ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        allRes = GenTable(sampleGOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = sum(score(resultFisher) <= 0.1))
        names(allRes) <- c("GO.ID", "Term", "Annotated", "Significant", "Expected", "classicFisher")
        ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        gg = genesInTerm(sampleGOdata)
        ggsig = lapply(gg, function(x) x[x %in% sigGenes])
        geneIds = lapply(ggsig, function(x) paste(x, collapse = ', '))
        gn = lapply(ggsig,  function(x) paste(mapvalues(x, from = Din[, 'ensembl_gene_id'], to = Din[, 'gene_name'], warn_missing = FALSE), collapse = ', '))
        gfc = lapply(ggsig, function(x) paste(mapvalues(x, from = Din[, 'ensembl_gene_id'], to = Din[, 'log2FoldChange'], warn_missing = FALSE), collapse = ', '))

        D_gogenes = data.frame(goid = names(gg), geneIds = unlist(geneIds), gn = unlist(gn), gfc = unlist(gfc))


        ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        D_out = merge(allRes, D_gogenes, by.x = 'GO.ID', by.y = 'goid', all.x = TRUE, sort = FALSE)
        D_out[, 'categoryName'] = GOTERM_categories[[go_cat]]
        D_out[, 'termName'] = paste(D_out[,'GO.ID'], D_out[, 'Term'], sep = '~')
        # D_out[,'percent'] = NA
        D_out[, 'listTotals'] = length(sigGenes)
        D_out[, 'listTotals_used'] = length(sigGenes(sampleGOdata))
        D_out[, 'popTotals'] = length(allGenes(sampleGOdata))
        D_out[, 'popTotals_used'] = numGenes(sampleGOdata)
        D_out[, 'foldEnrichment'] = (D_out[,'Significant']/D_out[,'listTotals_used']) / (D_out[, 'Annotated']/D_out[, 'popTotals_used'])

        # pull the column
        x <- D_out[, "classicFisher"]

        # coerce to character safely
        x_chr <- as.character(x)

        # parse: if it starts with "<", strip it and keep the number
        x_num <- suppressWarnings(as.numeric(sub("^\\\s*<\\\s*", "", x_chr)))

        # optional: clamp to a minimum positive value to avoid zeros/underflow issues
        x_num <- pmax(x_num, .Machine$double.xmin, na.rm = FALSE)

        D_out[, 'bonferroni'] = p.adjust(x_num, method = "bonferroni")
        D_out[, 'benjamini'] = p.adjust(x_num, method = "BY")
        D_out[, 'afdr'] = p.adjust(x_num, method = "fdr")

        ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        D_out = D_out[, c('categoryName', 'termName', 'Significant', 'classicFisher', 'geneIds',
                            'listTotals', 'listTotals_used', 'Annotated', 'popTotals', 'popTotals_used', 
                            'Expected',  'foldEnrichment', 'bonferroni', 'benjamini', 'afdr', "gn", "gfc")]

        names(D_out) <- c('categoryName', 'termName', 'listHits', 'classicFisher', 'geneIds',
                            'listTotals', 'listTotals_used', 'popHits', 'popTotals', 'popTotals_used',
                            'Expected', 'foldEnrichment', 'bonferroni', 'benjamini', 'afdr', 'genes name', 'log2fc')

        results[[GOTERM_categories[[go_cat]]]] = D_out
    }

    # save list as excel workbook
    write.xlsx(results, gsub("results.tsv","topGO.xlsx", fin), row.names = FALSE)
    # rbind and save as tsv
    results.tab = do.call(rbind, results[2:4])
    write.table(results.tab, gsub("results.tsv","topGO.tsv", fin), row.names = FALSE, sep = '\\t', quote = FALSE)
}
""",
    desc={
        "deseq2_output":"",
        "input_file":"",
    },
    container="mpgagebioinformatics/topgo:2.50.0",
    manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "4:00:00" }
)

cellplot=jawm.Process(
    name="cellplot",
    when=lambda p: ( not os.path.isfile( os.path.join( p.var["deseq2_output"], "annotated", p.var["inFile"].split( p.var["filetype"] )[0] )+p.var["category"]+".cellplot.pdf"  ) ) & \
                   ( not os.path.isfile( os.path.join( p.var["deseq2_output"], "annotated", p.var["inFile"].split( p.var["filetype"] )[0] )+p.var["category"]+".cellplot.touch" ) ) ,
    script="""\
#!/usr/bin/Rscript

# USAGE: `Rscript david_to_cellplot.R /beegfs/group_bit/data/projects/departments/Thomas_Langer/TL_Kai_RNAseq_5lines/adiff_output/Ka5R/wt_gfp.vs.cko_nlrp.DAVID.tsv tsv KEGG_PATHWAY 10`
# works with csv, tsv, txt and xlsx

# need to have CellPlot installed: devtools::install_github("dieterich-lab/CellPlot", build_vignettes = TRUE)
# and openxlsx: install.library('readxl')

rm(list = ls())

args<-commandArgs(TRUE)

#if (!require(devtools)) install.packages('devtools')
library(devtools)

#if (!require(readxl)) install.packages('readxl')
library(readxl)

#if (!require(CellPlot)) devtools::install_github("dieterich-lab/CellPlot")
library(CellPlot)

setwd("{{deseq2_output}}/annotated/")

# read in data, either excel or tsv
# reformat input data
inFile <- toString("{{inFile}}")
filetype <- toString("{{filetype}}")
category <- toString("{{category}}")
nterms <- as.numeric("{{nterms}}")

filetype_map <- c("xlsx" = 'xlsx',  'tsv' = '\t', 'csv' = ',', 'txt'=" ")

tmp_name=gsub(filetype, paste(category, '.cellplot.touch', sep = ''), inFile)


if(filetype == 'xlsx'){
    D <- read_excel(inFile, sheet = category)
    D <- as.data.frame(D)
  } else {
    D <- read.csv(inFile, header = TRUE, sep = filetype_map[filetype], as.is = TRUE)
}

D<-D[D["categoryName"] == category, ]
if ( nrow(D) == 0 ) {
  quit(save="no")
  file.create(tmp_name)
}

if( "classicFisher" %in%  names(D) ){
    x <- D[, "classicFisher"]

    # coerce to character safely
    x_chr <- as.character(x)

    # parse: if it starts with "<", strip it and keep the number
    x_num <- suppressWarnings(as.numeric(sub("^\\\s*<\\\s*", "", x_chr)))

    # optional: clamp to a minimum positive value to avoid zeros/underflow issues
    x_num <- pmax(x_num, .Machine$double.xmin, na.rm = FALSE)

  D$ease <- as.numeric(x_num)
}

D$ease <- as.numeric(as.character(D$ease))
D$foldEnrichment <- as.numeric(as.character(D$foldEnrichment))
D$listHits <- as.numeric(as.character(D$listHits))

# Added for handling cellplot NA values
# Handle log2fc properly
D$log2fc <- gsub("inf", "Inf", as.character(D$log2fc))

# if( !( "classicFisher" %in%  names(D) ) ){
# D$log2fc <- as.numeric(D$log2fc)
# }

print(1)
print( head(D) )

# Remove rows where `ease`, `foldEnrichment`, or `log2fc` are NA
D <- D[!is.na(D$ease) & !is.na(D$foldEnrichment) & !is.na(D$log2fc), ]

print(2)
print(head(D))

# If no valid data remains, exit
if (nrow(D) == 0) {
    message("No valid data available for plotting.")
    file.create(tmp_name)
    quit(save="no")
}
# End of addition for handling cellplot NA values

D <- D[order(D$ease),]

# subset to number of rows to plot..if specified number is larger than number of rows.
if (nterms >= nrow(D)) nterms = nrow(D)
D <- D[1:nterms,]

# log2FoldChange as list
D$log2fc <-lapply( gsub('inf', 'Inf', D$log2fc), function(x) as.numeric(as.character(unlist(strsplit(toString(x), ", ")))))

# cellplot    
x <- D

pdf(gsub(filetype, paste(category, '.cellplot.pdf', sep = ''), inFile))
cell.plot(x = setNames(-log10(D$ease), D$termName), 
              cells = D$log2fc, 
              main ="GO enrichment",
              xlab ="-log10(P.Value)", 
              x.mar = c(max(unlist(lapply(D$termName, function(x) nchar(x))))/100 + 0.1,0),
              key.n = 7, 
              y.mar = c(0.1, 0.1), 
              cex = 1.6, 
              cell.outer = 3, 
              bar.scale = .7, 
              space = .2)

dev.off()

# symplot
pdf(gsub(filetype, paste(category, '.symplot.pdf', sep = ''), inFile))
sym.plot(x = setNames(-log10(D$ease), D$termName), 
            cells = D$log2fc, 
            x.annotated = D$listHits, 
            main = "GO enrichment",
            key.lab = "-log10(P.Value)",
            x.mar = c(max(unlist(lapply(D$termName, function(x) nchar(x))))/100 + 0.1, 0),
            y.mar = c(0.2,0.1),
            key.n = 7, 
            cex = 1.6, 
            axis.cex = .8, 
            group.cex = .7) 

dev.off()
file.create(tmp_name)
""",
    desc={
        "deseq2_output":"",
        "inFile":"",
        "filetype":"",
        "category":"",
        "nterms":""
        
    },
    container="mpgagebioinformatics/cellplot:ea2dbc4",
    manager_slurm={ "-c": 8, "--mem": "16GB", "-t": "4:00:00" }
)

rcistarget=jawm.Process(
    name="rcistarget",
    when=lambda p: ( not os.path.isfile( os.path.join( p.var["deseq2_output"] , "annotated", p.var["input_file"].replace( ".results.tsv", ".RcisTarget.xlsx")  ) ) &  ( p.var["rcis_db"] != "" ) ),
    script="""\
#!/usr/bin/Rscript
library(RcisTarget)
library(openxlsx)
setwd("{{deseq2_output}}/annotated/")
# fout=stringr::str_split("{{input_file}}", ".results.tsv")[[1]][[1]]
# fout=paste0(fout,".RcisTarget.xlsx")
fout <- sub("\\\.results\\\.tsv$", ".RcisTarget.xlsx", "{{input_file}}")

# Load gene sets to analyze. e.g.:
deg = read.delim("{{input_file}}", as.is = TRUE)
sigGene = subset(deg, padj <= 0.05)[, 'gene_name']
sigGeneLists <- list(geneListName=sigGene)
motif_files = list(chip = c("{{rcis_db}}/chip_collection.feather", 
                            "{{rcis_db}}/chip_annotation.tsv"),
                  tf   = c("{{rcis_db}}/tf_collection.feather", 
                            "{{rcis_db}}/tf_annotation.tsv"))
results = list()
for(db in names(motif_files)){
    if(file.exists(motif_files[[db]][1])){
        ## 0. load motif data
        motifRankings <- importRankings(motif_files[[db]][1])
        motif_anno = read.delim(motif_files[[db]][2], as.is = TRUE)

        ## 1. Calculate AUC
        motifs_AUC <- calcAUC(sigGeneLists, motifRankings)

        ## 2. Select significant motifs, add TF annotation & format as table
        motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3)

        ## 3. Identify significant genes for each motif
        motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                          geneSets=sigGeneLists,
                                                          rankings=motifRankings)
        ## 4. Add motif name
        motifEnrichmentTable_wGenes = merge(motif_anno, motifEnrichmentTable_wGenes, by.x = 'X.motif_id', by.y = 'motif', all.y = TRUE)

        ## 5. Store data
        results[[db]] = motifEnrichmentTable_wGenes[order(motifEnrichmentTable_wGenes[, 'NES'], decreasing = TRUE), ]

        # make flattened table
        toFlat <- lapply(seq_len(nrow(motifEnrichmentTable_wGenes)), function(i) {
        genes_column <- unlist(strsplit(motifEnrichmentTable_wGenes[i, 'enrichedGenes'], ';', fixed = TRUE))
        data.frame(
            motifID   = motifEnrichmentTable_wGenes[i, 'X.motif_id'],
            motifName = motifEnrichmentTable_wGenes[i, 'gene_name'],
            targetGene= genes_column,
            NES       = motifEnrichmentTable_wGenes[i, 'NES'],
            stringsAsFactors = FALSE
        )
        })
        flat_df <- if (length(toFlat)) do.call(rbind, toFlat) else data.frame()
        results[[paste0('flat_', db)]] <- flat_df
    }
}
results <- lapply(results, function(x) as.data.frame(x, stringsAsFactors = FALSE))
wb <- createWorkbook()
for (nm in names(results)) {
  df <- as.data.frame(results[[nm]], stringsAsFactors = FALSE)
  addWorksheet(wb, nm)
  if (NROW(df) == 0 || NCOL(df) == 0) {
    writeData(wb, nm, data.frame(note = "(no data)"))
  } else {
    writeData(wb, nm, df)
    setColWidths(wb, nm, cols = 1:ncol(df), widths = "auto")
  }
}
saveWorkbook(wb, fout, overwrite = TRUE)
""",
    var={
        "rcis_db": ""
    },
    desc={
        "rcis_db": "Available from https://resources.aertslab.org/cistarget/",
        "deseq2_output":"",
        "input_file":""
    },
    container="mpgagebioinformatics/rcistarget:1.17.0",
    manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "4:00:00" }
)

# qc_plots=jawm.Process(
#     name="qc_plots",
#     when=lambda p: not os.path.isfile( os.path.join( p.var["deseq2_output"], "qc_plots" , "pca_all_samples.pdf"  ) ),
#     script="""\
# #!/bin/bash
# mkdir -p {{deseq2_output}}/qc_plots/
# QC_plots -o {{deseq2_output}}/qc_plots/ -de {{deseq2_output}}/ -s {{deseq2_output}}/samples_MasterTable.txt -t {{deseq2_output}} -sp {{spec}}
# """,
#     desc={
#         "deseq2_output":"",
#         "spec": ""
#     },
#     container="mpgagebioinformatics/rnaseq.python:3.8-8",
#     manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "4:00:00" }
# )


qc_plots=jawm.Process(
    name="qc_plots",
    when=lambda p: not os.path.isfile( os.path.join( p.var["deseq2_output"], "qc_plots" , "pca_all_samples.pdf"  ) ),
    script="""\
#!/usr/bin/env python3

from ast import arg
import os
import sys
import argparse

sys.stdout.flush()

# parser = argparse.ArgumentParser(description="Outputs multiple QC plots for differential expression data", \
# formatter_class = argparse.ArgumentDefaultsHelpFormatter)
# parser.add_argument("-o", "--output", help="/path/to/output/prefix")
# parser.add_argument("-de", "--DiffExpFolder", help="/path/to/diff_exp/folder")
# parser.add_argument("-s", "--sampleConditions", help="/path/to/sample/conditions/file")
# parser.add_argument("-t", "--topFolder", help="/path/to/top/folder/of/project")
# parser.add_argument("-sp", "--species", help="Supported values - c_albicans_sc5314,celegans,hsapiens,dmelanogaster,mmusculus,nfurzeri,scerevisiae")

# args = parser.parse_args()

import os
import AGEpy as age
import pandas as pd
import numpy as np
import matplotlib
import seaborn as sns

import matplotlib.pyplot as plt
import scipy
from  matplotlib.colors import LinearSegmentedColormap
from scipy import stats
from sklearn.cluster import KMeans
from sklearn import metrics
from scipy.spatial.distance import cdist
import multiprocessing as mp
from collections import defaultdict
from sklearn.metrics.pairwise import euclidean_distances
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import fcluster
import scipy.stats as stats
from matplotlib.pyplot import rc
from scipy.stats import hypergeom
import itertools
import sys
import statsmodels.stats.multitest as multi
from sklearn.decomposition import PCA
from sklearn import preprocessing
from itertools import cycle
import scipy.spatial.distance

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

#%matplotlib inline

###### SETTING VARIABLES #######

top = "{{deseq2_output}}/"
qcp = "{{deseq2_output}}/qc_plots/" 
deg = "{{deseq2_output}}/"
spec = "{{spec}}"
samples_treatment = "{{deseq2_output}}/samples_MasterTable.txt"

if not os.path.exists(top):
    raise ValueError("Error: Invalid path to Top Folder")

if not os.path.exists(deg):
    raise ValueError("Error: Invalid path to Differential Expression Folder")

if not os.path.exists(samples_treatment):
    raise ValueError("Error: Invalid path to Sample Conditions File")

if not os.path.isdir(qcp):
    os.makedirs( qcp )

##### PLOTTING #######
os.chdir(top)
os.getcwd()

df_original=pd.read_excel(deg+"/all_res_counts.xlsx", index_col=0, engine='openpyxl')
print(len(df_original))
df_original.head()

samp_df=pd.read_csv(samples_treatment, sep="\\t")
samp_df=samp_df.rename(columns={'Unnamed: 0':'Samples'})
group_col=samp_df.columns.tolist()[1]

samples_dict={}
conditions=list(set([s for s in samp_df[group_col].tolist()]))
for c in conditions:
    samples=samp_df.loc[samp_df[group_col] == c , 'Samples'].tolist()
    samples_dict[c]=samples


df_norm=np.log10(df_original+1)
df_norm.head()

font = {'family' : 'Serif',
            'weight' : 'semibold',
            'size'   : 12}
font_ = {'family' : 'Serif',
            'weight' : 'semibold',
            'size'   : 14}


#### GROUPED KDE ####

fig = plt.figure( frameon=False,figsize=(7,7))

for key in samples_dict.keys():
    dict_value=samples_dict[key]
    values=df_norm[dict_value].values.tolist()
    values=[item for sublist in values for item in sublist]
    values=[s for s in values if s != 0]
    sns.set()
    sns.kdeplot(values, shade=True, label=key)

plt.xlabel("log10(counts+1)",font_)
plt.ylabel("Density",font_)
plt.legend(loc='best', prop=font)
plt.title('Genes', font_)
plt.savefig(qcp+"/grouped.KDE.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### SAMPLE KDE ####

matplotlib.rc('font', **font)

fig = plt.figure( frameon=False,figsize=(7,7))

for f in df_norm.columns.tolist():
    values=[s for s in df_norm[f].tolist()]
    values=[s for s in values if s != 0]    
    sns.set()
    sns.kdeplot(values, shade=True, label=f)

plt.xlabel("log10(counts+1)",font_)
plt.ylabel("Density",font_)
plt.title('Genes', font_)
plt.legend(bbox_to_anchor=(1, 1), prop=font_)
plt.savefig(qcp+"/sample.KDE.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### GROUPED BAR PLOTS ####

matplotlib.rc('font', **font)

fig = plt.figure( frameon=False,figsize=(7,7))

medianprops = dict(color='black', linewidth=1)
#flierprops = dict(marker='D', markerfacecolor='black', markersize=6,linestyle='none')

cond_values=[]
for key in samples_dict.keys():
    dict_value=samples_dict[key]
    values=df_norm[dict_value].values.tolist()
    values=[item for sublist in values for item in sublist]
    values=[s for s in values if s != 0]
    cond_values.append(values)
  
sns.set()
bplot=plt.boxplot(cond_values, labels=conditions, patch_artist=True, showfliers=True, medianprops=medianprops)
plt.xticks(fontfamily='Serif')
plt.xlabel("condition",font_)
plt.ylabel("log10(counts+1)",font_)
palette = sns.color_palette("husl", len(conditions))
for patch, color in zip(bplot['boxes'], palette):
        patch.set_facecolor(color)

plt.savefig(qcp+"/grouped.barPlots.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### SAMPLE BAR PLOTS ####

samples=[col for col in df_norm.columns.tolist()]
values=[]

for col in df_norm.columns.tolist():
    tmp=df_norm[col].values.tolist()
    tmp=[s for s in tmp if s != 0]
    values.append(tmp)

matplotlib.rc('font', **font)

fig = plt.figure(frameon=False,figsize=(10,7))

sns.set()
bplot=plt.boxplot(values, labels=samples, patch_artist=True, showfliers=True, medianprops=medianprops)
plt.xticks(rotation=90, fontfamily='Serif')
plt.xlabel("condition",font_)
plt.ylabel("log10(counts+1)",font_)

palette = sns.color_palette("husl", len(samples))
for patch, color in zip(bplot['boxes'], palette):
        patch.set_facecolor(color)

plt.savefig(qcp+"/sample.barPlots.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')
plt.show()

#### SCATTER PLOT MATRIX ####

file_name=qcp+"/count.matrix.scatter.plot.pdf"

pdf = PdfPages(file_name)

w=len(samples_dict.keys())
l=len(samples_dict.keys())
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in samples_dict.keys():
    value_1=samples_dict[comp1]
    for comp2 in samples_dict.keys():
        value_2=samples_dict[comp2]
        if comp1 == comp2:            
            sns.set()
            ax = fig.add_subplot(l,w,i)
            values=df_norm[value_1].values.tolist()
            values=[item for sublist in values for item in sublist]
            values=[s for s in values if s != 0]
            sns.kdeplot(values)
            #ax.set_xlabel("log10(counts+1)",font_)
            ax.set_ylabel(comp1,font_)
            ax.set_title(comp2, font)
            i=i+1
        else:
            cols_y=value_1
            cols_x=value_2
            tmp = df_norm[value_1 + value_2]
            tmp = tmp[ tmp != 0]
            ax = fig.add_subplot(l,w,i)
            x=tmp[cols_x].mean(axis = 1)
            y=tmp[cols_y].mean(axis = 1)
            ax.scatter( x, y , color="tab:blue", alpha=0.5, s=1,)
            ax.set_xlabel("log10(counts+1)",font_)
            ax.set_ylabel(comp1+"\\n\\nlog10(counts+1)",font_)
            #ax.legend(loc='best', prop=font_)
            ax.set_title(comp2, font)
            i=i+1
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))
if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))
pdf.close()

#### DENDROGRAM ####

df_norm.T.index

plt.figure(figsize=(10, 10))
plt.ylabel('distance')
dendrogram(linkage(df_norm.T, 'ward'), labels=df_norm.T.index, leaf_rotation=90.0)
plt.savefig(qcp+"/sample.dendrogram.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

path=deg+"/annotated/"
data_files=os.listdir(path)
data_files=[s for s in data_files if ".results.tsv" in s]

conditions_new=[s.split("_",1)[1].split(".",1)[0].split("_vs_") for s in data_files]
conditions_new=[item for sublist in conditions_new for item in sublist]
conditions_new=list(set(conditions_new))
conditions_new

#### MA PLOTS ####

file_name=qcp+"/MA.plots.pdf"

pdf = PdfPages(file_name)

w=len(conditions_new)
l=len(conditions_new)
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            
            sns.set_style("whitegrid")
            ax = fig.add_subplot(l,w,i)
            
            ax.set_xlabel("log10(baseMean)",font_)
            ax.set_ylabel(comp1+"\\n\\nlog2FoldChange",font_)
            ax.set_title(comp2, font)
            
            i=i+1
            
        else:
            file=[s for s in data_files if (comp1+'.' in s or comp1+'_' in s) and (comp2+'.' in s or comp2+'_' in s) ]
            lfc_comp=[s.split(".results.tsv")[0] for s in file]
            df_tmp=pd.read_csv(path+file[0], sep="\\t")
            
            ax = fig.add_subplot(l,w,i)
            
            tmp=df_tmp[df_tmp["padj"]>=0.05]
            x=tmp["baseMean"].tolist()
            y=tmp["log2FoldChange"].tolist()
            x=[ np.log10(s+1) for s in x ]
            ax.scatter( x, y , color="k", alpha=0.5, s=1, label='NonSignificant')

            tmp=df_tmp[df_tmp["padj"]<0.05]
            x=tmp["baseMean"].tolist()
            y=tmp["log2FoldChange"].tolist()
            x=[ np.log10(s+1) for s in x ]
            ax.scatter( x, y , color="r", alpha=0.5, s=1, label='Significant')

            ax.set_xlabel("log10(baseMean)",font_)
            ax.set_ylabel(comp1+"\\n\\nlog2FoldChange("+lfc_comp[0]+")",font_)
            ax.legend(loc='best', prop=font_)
            ax.set_title(comp2, font)


            i=i+1
        
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))

if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))

pdf.close()

#### VOLCANO PLOTS ####

file_name=qcp+"/volcano.plots.pdf"

pdf = PdfPages(file_name)

w=len(conditions_new)
l=len(conditions_new)
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            
            sns.set_style("whitegrid")
            ax = fig.add_subplot(l,w,i)
            
            ax.set_xlabel("log2FoldChange",font_)
            ax.set_ylabel(comp1+"\\n\\n-log10(Adj.P.value)",font_)
            ax.set_title(comp2, font)
            
            i=i+1
            
        else:
            file=[s for s in data_files if (comp1+'.' in s or comp1+'_' in s) and (comp2+'.' in s or comp2+'_' in s) ]
            lfc_comp=[s.split(".results.tsv")[0] for s in file]
            df_tmp=pd.read_csv(path+file[0], sep="\\t")
            
            ax = fig.add_subplot(l,w,i)
            
            tmp=df_tmp[df_tmp["padj"]>=0.05]
            x=tmp["log2FoldChange"].tolist()
            y=tmp["padj"].tolist()
            y=[ -np.log10(s+1) for s in y ]
            ax.scatter( x, y , color="k", alpha=0.5, s=1, label='NonSignificant')

            tmp=df_tmp[df_tmp["padj"]<0.05]
            x=tmp["log2FoldChange"].tolist()
            y=tmp["padj"].tolist()
            y=[ -np.log10(s+1) for s in y ]
            ax.scatter( x, y , color="r", alpha=0.5, s=1, label='Significant')

            ax.set_xlabel("log2FoldChange("+lfc_comp[0]+")",font_)
            ax.set_ylabel(comp1+"\\n\\n-log10(Adj.P.value)",font_)
            ax.legend(loc='best', prop=font_)
            ax.set_title(comp2, font)

            i=i+1
        
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))

if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))

pdf.close()

#### P.VALUE DISTRIBUTION ####

file_name=qcp+"/p.value.dist.pdf"

pdf = PdfPages(file_name)

w=len(conditions_new)
l=len(conditions_new)
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            
            sns.set_style("whitegrid")
            ax = fig.add_subplot(l,w,i)
            
            ax.set_xlabel("P.Value",font_)
            ax.set_ylabel(comp1+"\\n\\nFrequency",font_)
            ax.set_title(comp2, font)
            
            i=i+1
            
        else:
            file=[s for s in data_files if (comp1+'.' in s or comp1+'_' in s) and (comp2+'.' in s or comp2+'_' in s) ]
            df_tmp=pd.read_csv(path+file[0], sep="\\t")
            
            ax = fig.add_subplot(l,w,i)
            
            df_tmp['pvalue'].hist(bins=30)
            
            ax.set_xlabel("P.Value",font_)
            ax.set_ylabel(comp1+"\\n\\nFrequency",font_)
            ax.set_title(comp2, font)

            i=i+1
        
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))

if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))

pdf.close()

#### Q.VALUE DISTRIBUTION ####

file_name=qcp+"/q.value.dist.pdf"

pdf = PdfPages(file_name)

w=len(conditions_new)
l=len(conditions_new)
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            
            sns.set_style("whitegrid")
            ax = fig.add_subplot(l,w,i)
            
            ax.set_xlabel("Q.Value",font_)
            ax.set_ylabel(comp1+"\\n\\nFrequency",font_)
            ax.set_title(comp2, font)
            
            i=i+1
            
        else:
            file=[s for s in data_files if (comp1+'.' in s or comp1+'_' in s) and (comp2+'.' in s or comp2+'_' in s) ]
            df_tmp=pd.read_csv(path+file[0], sep="\\t")
            
            ax = fig.add_subplot(l,w,i)
            
            df_tmp['padj'].hist(bins=30)
            
            ax.set_xlabel("Q.Value",font_)
            ax.set_ylabel(comp1+"\\n\\nFrequency",font_)
            ax.set_title(comp2, font)

            i=i+1
        
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))

if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))

pdf.close()

#### GENESET LEVEL PLOTS ####

genes=[]
for f in data_files:
    tmp=pd.read_csv(path+f, sep="\\t")
    tmp_genes=tmp.loc[tmp['padj'] < 0.05 , 'ensembl_gene_id'].tolist()
    for g in tmp_genes:
        if g not in genes:
            genes.append(g)
            
len(genes)

df_heat=df_norm.loc[df_norm.index.isin(genes),]
df_heat.head()

#print(conditions)

for key in samples_dict.keys():
    value=samples_dict[key]
    df_heat['mean('+key+')']=df_heat[value].mean(axis=1)
df_heat.head()

#### GROUPED HEATMAP ####

matplotlib.rc('font', **font)

#fig = plt.figure( frameon=False,figsize=(14,14))
cm=sns.clustermap(df_heat[[s for s in df_heat.columns.tolist() if 'mean' in s]],figsize=(14,14))
cm.ax_row_dendrogram.set_visible(False)
cm.ax_col_dendrogram.set_visible(False)
plt.savefig(qcp+"/grouped.heatMap.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### SAMPLE HEATMAP ####

matplotlib.rc('font', **font)

#fig = plt.figure( frameon=False,figsize=(14,14))
cm=sns.clustermap(df_heat[[s for s in df_heat.columns.tolist() if 'mean' not in s]],figsize=(14,16))
cm.ax_row_dendrogram.set_visible(False)
cm.ax_col_dendrogram.set_visible(False)
plt.savefig(qcp+"/sample.heatMap.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### SIGNIFICANT FEATURES MATRIX ####

df_sig=pd.DataFrame(columns=conditions_new, index=conditions_new)

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            sig_fearures=0
            df_sig.loc[comp1,comp2]=sig_fearures
        else:
            file=[s for s in data_files if (comp1+'.' in s or comp1+'_' in s) and (comp2+'.' in s or comp2+'_' in s) ]
            tmp=pd.read_csv(path+file[0], sep="\\t")
            sig_fearures=len(tmp.loc[tmp['padj'] < 0.05,])
            df_sig.loc[comp1,comp2]=sig_fearures
        
# df_sig

mask = np.triu(np.ones_like(df_sig.astype(float), dtype=np.bool))

matplotlib.rc('font', **font)

fig = plt.figure( frameon=False,figsize=(9,7))

sns.heatmap(df_sig.astype(float),annot=True,fmt='.10g', mask=mask)

plt.title('No of significant features', font_)
plt.yticks(rotation=0) 

plt.savefig(qcp+"/sigFeatures.matrix.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

df_corr=df_norm.copy()
for key in samples_dict.keys():
    value=samples_dict[key]
    df_corr['mean('+key+')']=df_corr[value].mean(axis=1)

df_corr.head()

#### GROUP DISTANCE MATRIX ####

dist=metrics.pairwise_distances(np.array(df_corr[[s for s in df_corr.columns.tolist() if 'mean' in s]].T))

df_dist=pd.DataFrame(dist, columns=[s.split("(")[1].strip(")") for s in df_corr.columns.tolist() if 'mean' in s], index=[s.split("(")[1].strip(")") for s in df_corr.columns.tolist() if 'mean' in s])
matplotlib.rc('font', **font)

annot = {'family' : 'Serif',
         'weight' : 'semibold',
         'size'   : 12}

fig = plt.figure( frameon=False,figsize=(12,10))

sns.heatmap(df_dist.T,annot=True,fmt='.5g',annot_kws=annot)

plt.title('Group Distance Matrix', font_)

plt.savefig(qcp+"/group.distance.matrix.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')   
plt.show()

#### SAMPLE DISTANCE MATRIX ####

dist_all=metrics.pairwise_distances(np.array(df_corr[[s for s in df_corr.columns.tolist() if 'mean' not in s]].T))

df_dist_all=pd.DataFrame(dist_all, columns=[s for s in df_corr.columns.tolist() if 'mean' not in s], index=[s for s in df_corr.columns.tolist() if 'mean' not in s])
matplotlib.rc('font', **font)

annot = {'family' : 'Serif',
         'weight' : 'semibold',
         'size'   : 10}

fig = plt.figure( frameon=False,figsize=(12,10))

sns.heatmap(df_dist_all.T,annot=True,fmt='.3g',annot_kws=annot)

plt.title('Sample Distance Matrix', font_)

plt.savefig(qcp+"/sample.distance.matrix.pdf", dpi=600,bbox_inches='tight', pad_inches=0.1,format='pdf')
plt.show()

#### PCA ####

file_name=qcp+"/pca.pdf"

pdf = PdfPages(file_name)

w=len(conditions_new)
l=len(conditions_new)
i=1

fig = plt.figure(figsize=(w*12,l*12))

for comp1 in conditions_new:
    for comp2 in conditions_new:
        if comp1 == comp2:
            
            sns.set_style("darkgrid")
            ax = fig.add_subplot(l,w,i)
            
            ax.set_xlabel("Component1",font_)
            ax.set_ylabel(comp1+"\\n\\nComponent2",font_)
            ax.set_title(comp2, font)
            
            i=i+1
            
        else:
            file=[s for s in data_files if (comp1+'.' in s or comp1+'_' in s) and (comp2+'.' in s or comp2+'_' in s) ]
            df_tmp=pd.read_csv(path+file[0], sep="\t")
            
            ax = fig.add_subplot(l,w,i)
            
            data=df_tmp.copy()
            
            #cols=[s for s in data.columns if comp1 in s or comp2 in s]
            #data_pca=data[cols]
            if spec != "c_albicans_sc5314":
                cols_to_exclude=['ensembl_gene_id','gene_name', 'baseMean','log2FoldChange','lfcSE','pvalue','padj','gene_biotype','GO_id','GO_term']
            else:
                cols_to_exclude=['ensembl_gene_id','gene_name', 'baseMean','log2FoldChange','lfcSE','pvalue','padj']
            print(data)
            data_pca=data.drop(cols_to_exclude, axis=1)
            data_pca=np.log10(data_pca + 1)
            data_pca=data_pca.set_index(data['ensembl_gene_id'])
            #data_pca.head()

            df_pca=data_pca.T.reset_index()
            #df_pca

            #print(df_pca.shape)

            df_pca.set_index('index', inplace=True)
            df_pca.index.names = ['Sample']

            pca = PCA(copy=True, iterated_power='auto', n_components=2, random_state=None,svd_solver='auto', tol=0.0, whiten=False)

            #scaling the values
            df_pca_scaled = preprocessing.scale(df_pca, axis = 1)

            projected=pca.fit_transform(df_pca_scaled)

            #print(pca.explained_variance_ratio_)

            tmp=pd.DataFrame(projected)
            tmp.rename(columns={0: 'Component 1',1: 'Component 2'}, inplace=True)
            tmp.to_excel(qcp+"/"+file[0].replace(".results.tsv","_pca.xlsx"), index=None)
            #tmp.head()

            final_pca=pd.merge(df_pca.reset_index(),tmp,left_index=True,right_index=True)
            #final_pca.head()

            for s in final_pca['Sample'].tolist():
                for key in samples_dict.keys():
                    value=samples_dict[key]
                    if s in value:
                        final_pca.loc[final_pca['Sample'] == s,'expt'] = key

            font = {'family' : 'serif',
                    'weight' : 'bold',
                    'size'   : 12}


            color_gen = cycle(('blue', 'red', 'green', 'pink','yellow'))
            
            for lab in set(final_pca['expt']):
                plt.scatter(final_pca.loc[final_pca['expt'] == lab, 'Component 1'], 
                            final_pca.loc[final_pca['expt'] == lab, 'Component 2'], 
                            c=next(color_gen),
                            label=lab)

            ax.set_xlabel('Component 1  - ' + str(pca.explained_variance_ratio_[0]*100)+ " % ", fontdict=font)
            ax.set_ylabel(comp1+"\\n\\nComponent 2  - "+ str(pca.explained_variance_ratio_[1]*100)+ "  % ", fontdict=font)
            #ax.set_xticks(fontsize=14)
            #ax.set_yticks(fontsize=14)
            ax.legend(loc='best', fontsize=14, prop=font_)
            ax.set_title(comp2, fontdict=font)

            i=i+1
        
        if i == w*l+1:
            plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
	    #plt.show()
            #plt.close()
            i=1
            fig = plt.figure(figsize=(w*12,l*12))

if i != 1:
    plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')
    #plt.show()
    #plt.close()
    i=1
    fig = plt.figure(figsize=(w*12,l*12))

pdf.close()



### PCA all ####

file_name=qcp+"/pca_all_samples.pdf"

pdf = PdfPages(file_name)

plt.figure(figsize=(14,14))

#sns.set_style("darkgrid")

pca_data = df_heat[[s for s in df_heat.columns.tolist() if 'mean' not in s]]
df_pca=pca_data.T.reset_index()

df_pca.set_index('index', inplace=True)
df_pca.index.names = ['Sample']

pca = PCA(copy=True, iterated_power='auto', n_components=2, random_state=None,svd_solver='auto', tol=0.0, whiten=False)

#scaling the values
df_pca_scaled = preprocessing.scale(df_pca, axis = 1)

projected=pca.fit_transform(df_pca_scaled)
#print(pca.explained_variance_ratio_)

tmp=pd.DataFrame(projected)
tmp.rename(columns={0: 'Component 1',1: 'Component 2'}, inplace=True)
#tmp.to_excel(qcp+"/pca_comp_all_samples.xlsx", index=True, header=True)
#tmp.head()

final_pca=pd.merge(df_pca.reset_index(),tmp,left_index=True,right_index=True)
#final_pca.head()

for s in final_pca['Sample'].tolist():
    for key in samples_dict.keys():
        value=samples_dict[key]
        if s in value:
            final_pca.loc[final_pca['Sample'] == s,'expt'] = key

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 12}

color_gen = cycle(('blue', 'red', 'green', 'pink','yellow'))
color_gen = sns.color_palette(None, len(conditions))

plot_data = final_pca[["Sample", "Component 1", "Component 2", "expt"]]
plot_data.to_excel(qcp+"/pca_comp_all_samples.xlsx", index=None)

# map group color
plot_data['color'] = plot_data['expt'].map(dict(zip(set(plot_data["expt"]), color_gen)))

for lab in set(plot_data['expt']):
    scatter = plt.scatter(plot_data.loc[plot_data['expt'] == lab, 'Component 1'],
    plot_data.loc[plot_data['expt'] == lab, 'Component 2'],
    c=plot_data.loc[plot_data['expt'] == lab, 'color'],
    label=lab)

plt.legend(handles=scatter.legend_elements()[0], labels=set(plot_data["expt"]), title = 'Group')

plt.xlabel('Component 1  - ' + str(pca.explained_variance_ratio_[0]*100)+ " % ", fontdict=font)
plt.ylabel('Component 2  - ' + str(pca.explained_variance_ratio_[1]*100)+ "  % ", fontdict=font)
plt.title("PCA of all Samples", font)

plt.savefig(pdf, dpi=300, bbox_inches='tight', pad_inches=0.1,format='pdf')

pdf.close()
""",
    desc={
        "deseq2_output":"",
        "spec": ""
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "4:00:00" }
)



get_ip=jawm.Process(
    name="get_ip",
    when=lambda p: ( not os.path.isfile( os.path.join( p.var["deseq2_output"] , "annotated", "string.done"  ) ) &  ( p.var["cytoscape_host"] != "" ) ) ,
    script="""\
#!/bin/bash
while [[ ! -f {{cytoscape_host}} ]] ; do 
  echo "waiting for cytoscape to be available"
  sleep 3$((RANDOM % 9))
done
mv {{cytoscape_host}} {{cytoscape_host}}_inuse
touch {{deseq2_output}}/annotated/string.running
""",
    desc={
        "cytoscape_host":""
    },
)

# using string over cytsocape is deprecated 
# a new module without cytoscape is currently under development

string=jawm.Process(
    name="string",
    when=lambda p: ( not os.path.isfile( os.path.join( p.var["deseq2_output"] , "annotated", "string.done"  ) ) &  ( p.var["cytoscape_host"] != "" ) ),
    script="""\
#!/usr/local/bin/python
import pandas as pd
import numpy as np
import AGEpy as age
import sys
import os
from pathlib import Path
from py2cytoscape import cyrest
from py2cytoscape.cyrest.base import *
import paramiko
from time import sleep
import matplotlib
import matplotlib.pyplot as plt
import tempfile
import traceback

try:

    ################# in values ################################
    with open("{{cytoscape_host}}_inuse" , "r") as hostfile:
        host=hostfile.readlines()[0].split("\\n")[0]
    print(f"host::{host}::")
    species="{{species}}"
    biomarthost="{{biomart_host}}"
    cytoscape=cyrest.cyclient(host=host)
    cytoscape.version()
    cytoscape.command.run("/tmp/cyscript")

    ###########################################################
    taxons={"caenorhabditis elegans":"6239","drosophila melanogaster":"7227",\
            "mus musculus":"10090","homo sapiens":"9606", "saccharomyces cerevisiae": "4932", "nothobranchius furzeri": "105023"}
    tags={"caenorhabditis elegans":"CEL","drosophila melanogaster":"DMEL",\
            "mus musculus":"MUS","homo sapiens":"HSA"}
    taxon_id=taxons[species]
    ###########################################################
    aging_genes = []
    ### ATTENTION ### if you are using yeast, you will need to uncomment the follwing lines
    input_files=os.listdir("{{deseq2_output}}/annotated/")
    input_files=[s for s in input_files if ".results.tsv" in s ]
    input_files=[ os.path.join("{{deseq2_output}}/annotated/",s) for s in input_files if ".results.tsv" in s ]
    python_output="/".join(input_files[0].split("/")[:-1])
    if species in tags.keys():
        organismtag=tags[species]

        if not os.path.isfile(python_output+"/homdf.txt"):
            print("Could not find ageing evidence table. Using biomart to create one.")
            sys.stdout.flush()
            homdf,HSA,MUS,CEL,DMEL=age.FilterGOstring(host=biomarthost)
            homdf.to_csv(python_output+"/homdf.txt", index=None,sep="\t")
        else:
            print("Found existing ageing evidence table.")
            sys.stdout.flush()
            homdf=pd.read_csv(python_output+"/homdf.txt", sep="\t")
            aging_genes=homdf[[organismtag+"_ensembl_gene_id","evidence"]].dropna()
            aging_genes=aging_genes[aging_genes[organismtag+"_ensembl_gene_id"]!="None"]
            aging_genes=aging_genes[organismtag+"_ensembl_gene_id"].tolist()
    ### till here
    ###########################################################
    input_files=os.listdir("{{deseq2_output}}/annotated/")
    input_files=[s for s in input_files if ".results.tsv" in s ]
    input_files=[ os.path.join("{{deseq2_output}}/annotated/",s) for s in input_files if ".results.tsv" in s ]
    for fin in input_files:
        python_output="/".join(fin.split("/")[:-1])
        target=fin.replace("results.tsv","cytoscape")
        if os.path.isfile(target+".cys"):
            continue


        dfin=pd.read_csv(fin, sep="\\t")
        cytoscape=cyrest.cyclient(host=host)
        cytoscape.version()
        cytoscape.session.new()
        # cytoscape.vizmap.apply(styles="default")
        # Annotate aging evindence
        def CheckEvidence(x,aging_genes=aging_genes):
            if x in aging_genes:
                res="aging_gene"
            else:
                res="no"
            return res
        ### also comment this line
        dfin["evidence"]=dfin["ensembl_gene_id"].apply(lambda x:CheckEvidence(x) )
        dfin["baseMean"]=dfin["baseMean"].apply(lambda x: np.log10(x))
        qdf=dfin[dfin["padj"]<0.05]
        if qdf.shape[0] == 0:
            sys.exit()
        qdf=qdf.sort_values(by=["padj"],ascending=True)
        query_genes=qdf["ensembl_gene_id"].tolist()[:1000]
        limit=int(len(query_genes)*.25)
        response=api("string", "protein query",\
                                {"query":",".join(query_genes),\
                                "cutoff":str(0.4),\
                                "species":species,\
                                "limit":str(limit),\
                                "taxonID":taxon_id},\
                                host=host, port="1234")
        cytoscape.layout.force_directed(defaultSpringCoefficient=".000004", defaultSpringLength="5")
        defaults_dic={"NODE_SHAPE":"ellipse",\
                        "NODE_SIZE":"60",\
                        "NODE_FILL_COLOR":"#AAAAAA",\
                        "EDGE_TRANSPARENCY":"120"}
        defaults_list=cytoscape.vizmap.simple_defaults(defaults_dic)
        NODE_LABEL=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_LABEL",\
                                                        mappingType="passthrough",\
                                                        mappingColumn="display name")
        cytoscape.vizmap.create_style(title="dataStyle",\
                                        defaults=defaults_list,\
                                        mappings=[NODE_LABEL])
        sleep(4)
        cytoscape.vizmap.apply(styles="dataStyle")
        uploadtable=dfin[dfin["padj"]<0.05][["ensembl_gene_id","baseMean","log2FoldChange","evidence"]].dropna()
        # uploadtable=dfin[dfin["padj"]<0.05][["ensembl_gene_id","baseMean","log2FoldChange"]].dropna() ### use this line if you are using yeast
        cytoscape.table.loadTableData(uploadtable,df_key="ensembl_gene_id",table_key_column="query term")
        sleep(10)
        cmap = matplotlib.cm.get_cmap("bwr")
        norm = matplotlib.colors.Normalize(vmin=-4, vmax=4)
        min_color=matplotlib.colors.rgb2hex(cmap(norm(-4)))
        center_color=matplotlib.colors.rgb2hex(cmap(norm(0)))
        max_color=matplotlib.colors.rgb2hex(cmap(norm(4)))  
        NODE_FILL_COLOR=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_FILL_COLOR",mappingType="continuous",\
                                                            mappingColumn="log2FoldChange",\
                                                        lower=[-4,min_color],\
                                                            center=[0.0,center_color],\
                                                            upper=[4,max_color])
        ### do not do this if you are using yeast ...
        # apply diamond shape and increase node size to nodes with aging evidence
        NODE_SHAPE=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_SHAPE",mappingType="discrete",mappingColumn="evidence",\
                                                        discrete=[ ["aging_gene","no"], ["DIAMOND", "ellipse"] ])
        NODE_SIZE=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_SIZE",mappingType="discrete",mappingColumn="evidence",\
                                                    discrete=[ ["aging_gene","no"], ["100.0","60.0"] ])
        ###
        cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_SIZE,NODE_SHAPE,NODE_FILL_COLOR])
        # cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_FILL_COLOR]) # if using yeast
        cytoscape.vizmap.apply(styles="dataStyle")
        network = "current"
        namespace='default'
        PARAMS=set_param(["columnList","namespace","network"],["SUID",namespace,network])
        network=api(namespace="network", command="get attribute",PARAMS=PARAMS, host=host,port='1234',version='v1')
        network=int(network[0]["SUID"])
        basemean = cytoscape.table.getTable(table="node",columns=["baseMean"], network = network)
        min_NormInt = min(basemean.dropna()["baseMean"].tolist())
        max_NormInt = max(basemean.dropna()["baseMean"].tolist())
        cent_NormInt = np.mean([min_NormInt,max_NormInt])
        cmap = matplotlib.cm.get_cmap("Reds")
        norm = matplotlib.colors.Normalize(vmin=min_NormInt, vmax=max_NormInt)
        min_color=matplotlib.colors.rgb2hex(cmap(norm(np.mean([min_NormInt,max_NormInt]))))
        center_color=matplotlib.colors.rgb2hex(cmap(norm(cent_NormInt)))
        max_color=matplotlib.colors.rgb2hex(cmap(norm(max_NormInt)))  
        NODE_BORDER_PAINT=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_BORDER_PAINT",\
                                                            mappingType="continuous",\
                                                            mappingColumn="baseMean",\
                                                            lower=[min_NormInt,min_color],\
                                                            center=[np.mean([min_NormInt,max_NormInt]),center_color],\
                                                            upper=[max_NormInt,max_color])
        cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_BORDER_PAINT])
        NODE_BORDER_WIDTH=cytoscape.vizmap.mapVisualProperty(visualProperty="NODE_BORDER_WIDTH",\
                                                            mappingType="continuous",\
                                                            mappingColumn="baseMean",\
                                                            lower=[min_NormInt,2],\
                                                            center=[np.mean([min_NormInt,max_NormInt]),4],\
                                                            upper=[max_NormInt,8])
        cytoscape.vizmap.update_style("dataStyle",mappings=[NODE_BORDER_WIDTH])
        cytoscape.vizmap.apply(styles="dataStyle")
        cytoscape.network.rename(name="main String network")
        cytoscape.network.select(edgeList="all", extendEdges="true")
        cytoscape.network.create(source="current",nodeList="selected")
        cytoscape.network.rename(name="main String network (edges only)")
        cytoscape.network.set_current(network="main String network (edges only)")
        log2FoldChange = cytoscape.table.getTable(table="node",columns=["log2FoldChange"])
        if int(len(log2FoldChange)*.10) > 0:
            log2FoldChange["log2FoldChange"]=log2FoldChange["log2FoldChange"].apply(lambda x: abs(x))
            log2FoldChange=log2FoldChange.sort_values(by=["log2FoldChange"],ascending=False)
            top_nodes=log2FoldChange.index.tolist()[:int(len(log2FoldChange)*.10)]
            cytoscape.network.set_current(network="main String network (edges only)")
            cytoscape.network.select(nodeList="name:"+",".join(top_nodes))
            cytoscape.network.select(firstNeighbors="any",network="current")
            sleep(5)
            cytoscape.network.create(source="current",nodeList="selected")
            cytoscape.network.rename(name="top "+str(int(len(log2FoldChange)*.10))+" changed firstNeighbors")
        def MAKETMP():
            (fd, f) = tempfile.mkstemp()
            f="/tmp/"+f.split("/")[-1]
            return f
        cys=MAKETMP()
        cyjs=MAKETMP()
        main_png=MAKETMP()
        main_pdf=MAKETMP()
        edg_png=MAKETMP()
        edg_pdf=MAKETMP()
        neig_png=MAKETMP()
        neig_pdf=MAKETMP()
        cytoscape.session.save_as(session_file=cys)
        cytoscape.network.export(options="CYJS",OutputFile=cyjs)
        cytoscape.network.set_current(network="main String network")
        cytoscape.network.deselect(edgeList="all",nodeList="all")
        cytoscape.view.export(options="PNG",outputFile=main_png)
        cytoscape.view.export(options="PDF",outputFile=main_pdf)
        cytoscape.network.set_current(network="main String network (edges only)")
        cytoscape.network.deselect(edgeList="all",nodeList="all")
        cytoscape.view.export(options="PNG",outputFile=edg_png)
        cytoscape.view.export(options="PDF",outputFile=edg_pdf)
        if int(len(log2FoldChange)*.10) > 0:
            cytoscape.network.set_current(network="top "+str(int(len(log2FoldChange)*.10))+" changed firstNeighbors")
            cytoscape.network.deselect(edgeList="all",nodeList="all")
            sleep(5)
            cytoscape.view.export(options="PNG",outputFile=neig_png)
            cytoscape.view.export(options="PDF",outputFile=neig_pdf)
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(host, username="flaski")
        ftp_client=ssh.open_sftp()
        for f, extension, local in zip([cys,cyjs,main_png,main_pdf,edg_png,edg_pdf,neig_png,neig_pdf],\
                                        [".cys",".cyjs",".png",".pdf",".png",".pdf",".png",".pdf" ],\
                                        [target+".cys",target+".cyjs",target+".main.png",target+".main.pdf",\
                                        target+".main.edges.png",target+".main.edges.pdf",\
                                        target+".topFirstNeighbors.png",target+".topFirstNeighbors.pdf"]):
            try:
                ftp_client.get(f+extension,local)
                ssh_stdin, ssh_stdout, ssh_stderr = ssh.exec_command("rm "+f+extension )
            except:
                print("No "+local)
                sys.stdout.flush()
        print(f"Done with cytoscape for {fin}.")
except Exception as e:
    os.rename("{{cytoscape_host}}_inuse", "{{cytoscape_host}}")
    print("Error:", e)
    traceback.print_exc()
    exit(1)

file_path = Path("{{deseq2_output}}/annotated/string.done" )
file_path.touch()
""",
    desc={
        "cytoscape_host":"Path to a txt file containing one line with the IP address of your running cytoscape instance.",
        "biomart_host": "",
        "species": "",
        "deseq2_output":""
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "20GB", "-t": "4:00:00" }
)

release_ip=jawm.Process(
    name="release_ip",
    when=lambda p: ( os.path.isfile( os.path.join( p.var["deseq2_output"] , "annotated", "string.running"  ) ) &  ( p.var["cytoscape_host"] != "" ) ),
    script="""\
#!/bin/bash
if [[ -f {{cytoscape_host}}_inuse ]] ; then mv {{cytoscape_host}}_inuse {{cytoscape_host}} ; fi
rm -rf {{deseq2_output}}/annotated/string.running
""",
    desc={
        "cytoscape_host": "",
        "deseq2_output":""
    },
    container=""
)


def report_files(deseq2_output) :
    report_paths={}
    dic={ 
        deseq2_output :{
            "deseq2": ["all.samples.raw.counts.tsv" ]
        },
        os.path.join( deseq2_output, "annotated") : { 
            "deseq2":[ "*.results.xlsx", "masterTable_annotated.xlsx", "significant.xlsx", "tpm.xlsx", "*vsd.counts.xlsx"] , 
            "david":[ "*DAVID*xlsx", "*DAVID*cellplot*pdf" ], 
            "rcistarget": "*.RcisTarget.*",
            "topgo":[ "*.topGO.xlsx", "*.topGO.*.cellplot.pdf" ] 
            },
        os.path.join( deseq2_output, "qc_plots" ) : { 
            "qc_plots":[ "*.pdf", "*.xlsx" ],
            }
        }

    for path in dic :
        directory = Path( path )
        for folder in dic[path] :
            pattern = dic[path][folder]
            if isinstance( pattern, str ):
                pattern=[ pattern ]
            files=[  ]
            for p in pattern :
                files=files+[ f.resolve() for f in directory.glob( p ) ]
            # print(str(path),"; ", str(folder),"; ", dic[path][folder] )
            if files :
                # print("\t",",".join( [ str(f) for f in files ]))
                if folder not in report_paths :
                    report_paths[folder]=files
                else:
                    report_paths[folder]=report_paths[folder]+files

    return report_paths


if __name__ == "__main__":
    import sys
    from jawm.utils import workflow
    import glob

    workflows, var, args, unknown_args = jawm.utils.parse_arguments(['main','deseq2','test', "cistarget"])

    if workflow(["main","deseq2","test","cistarget"], workflows):

        tx2gene.execute( )

        if not workflow("test", workflows):
            # we can not run the biomart on github
            annotations.execute()

        parse_submission.execute()
        tx2gene_proc.execute(tx2gene.hash)

        jawm.Process.wait([ tx2gene_proc.hash, parse_submission.hash ])

        base = deseq2.var["deseq2_output"]
        tests = glob.glob(os.path.join(base, "*input.tsv"))

        deseq2_jobs=[]

        for file in tests:

            # clone the processes required for each file
            deseq2_=deseq2.clone()

            deseq2_.var["input_file"]=os.path.basename(file)
            deseq2_.execute()

            deseq2_jobs.append( deseq2_.hash )
        
        mastertable.execute( deseq2_jobs )

        ########### TEST
        if workflow("test", workflows):
            import pandas as pd
            # for the test workflow we might also do something more
            with open( os.path.join(var["deseq2_output"], "all_results_stats.tsv"), 'r') as out:
                infile = os.path.join(var["deseq2_output"], "all_results_stats.tsv")
                outfile = os.path.join(var["deseq2_output"], "all_results_stats.rounded.tsv")
                df = pd.read_csv(infile, sep="\t")
                # Round all float columns to 6 significant digits
                float_cols = df.select_dtypes(include="float").columns
                df[float_cols] = df[float_cols].map(lambda x: float(f"{x:.6f}") if pd.notna(x) else x)
                df.to_csv(outfile, sep="\t", index=False, float_format="%.6f", na_rep="NA")

            with open( os.path.join(var["deseq2_output"], "all_results_stats.rounded.tsv"), 'r') as out:
                print( "".join(out.readlines()[:2] + out.readlines()[-1:])  )

            # we can not run the remaining part of the workflow on github
            # so we stop it here
            with open( os.path.join(var["deseq2_output"], "test.txt"), 'w') as out:
                out.write("Test completed.")
            print("Test completed.")
            sys.exit(0)
        ########### TEST 

        annotator.execute( mastertable.hash )

        # wait for annotator to complete before listing input files and starting david
        # as they are generated by annotator
        jawm.Process.wait( annotator.hash )

        david.execute( )

        base = os.path.join(deseq2.var["deseq2_output"], "annotated")
        tests_results = glob.glob(os.path.join(base, "*results.tsv"))
        topgo_jobs=[]
        for file in tests_results:
            topgo_=topgo.clone()
            topgo_.var["input_file"]=os.path.basename(file)
            topgo_.execute( )
            topgo_jobs.append( topgo_.hash )

        jawm.Process.wait( topgo_jobs + [ david.hash ] )

        base = os.path.join(deseq2.var["deseq2_output"],  "annotated")
        david_files=glob.glob(os.path.join(base, "*.DAVID.tsv"))
        togo_files=glob.glob(os.path.join(base, "*.topGO.tsv"))

        cellplot_jobs=[]
        for file in david_files:
            cellplot_=cellplot.clone()

            cellplot_.var["map.inFile"]=file
            cellplot_.var["category"]="GOTERM_BP_FAT"

            cellplot__=cellplot_.clone()
            cellplot__.var["category"]="KEGG_PATHWAY"

            cellplot_.execute()
            cellplot__.execute()

            cellplot_jobs.append( cellplot_.hash )
            cellplot_jobs.append( cellplot__.hash )

        for file in togo_files:

            cellplot_=cellplot.clone()

            cellplot_.var["map.inFile"]=file
            cellplot_.var["category"]="GOTERM_BP"
            cellplot_.execute()

            cellplot_jobs.append( cellplot_.hash )

        qc_plots.execute( )

        # deprecated / see above
        # get_ip.execute( qc_plots.hash )
        # string.execute( get_ip.hash )
        # release_ip.execute( string.hash )
        

    if workflow("cistarget", workflows):

        rcistarget_jobs=[]
        for file in tests_results:
            rcistarget_=rcistarget.clone()
            rcistarget_.var["input_file"]=os.path.basename(file)
            rcistarget_.execute( )
            rcistarget_jobs.append( rcistarget_.hash )

    jawm.Process.wait()


