import jawm
import os

annotations=jawm.Process(
    name="annotations",
    when=lambda p: not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" ,"annotated",  "biotypes_go.txt") ) ,
    script="""\
#!/usr/local/bin/python
import pandas as pd
import os
import AGEpy as age
import shutil

if not os.path.isdir("{{project_folder}}/deseq2_output/annotated/") :
    os.makedirs("{{project_folder}}/deseq2_output/annotated/")

if "{{biomart_host}}" != "None":
    if not os.path.isfile("{{biotypes_go}}"):
        from biomart import BiomartServer
        attributes=["ensembl_gene_id","external_gene_name","go_id","name_1006"]
        try:
            server = BiomartServer( "{{biomart_host}}" )
            organism=server.datasets["{{biomart_dataset}}"]
            response=organism.search({"attributes":attributes})
            response=response.content.decode().split("\\n")
            response=[s.split("\\t") for s in response ]
            bio_go=pd.DataFrame(response,columns=attributes)
            bio_go=bio_go.sort_values(by=attributes, ascending=True )
            bio_go.to_csv("{{biotypes_go}}".replace("biotypes_go.txt", "biotypes_go_raw_topgo.txt"), index=None, sep="\\t")
            bio_go.to_csv("{{project_folder}}/deseq2_output/annotated/biotypes_go_raw_topgo.txt", index=None, sep="\\t")
            bio_go=bio_go[["ensembl_gene_id","go_id","name_1006"]]
            bio_go.to_csv("{{biotypes_go}}".replace("biotypes_go.txt", "biotypes_go_raw.txt"), index=None, sep="\\t")
            bio_go.to_csv("{{project_folder}}/deseq2_output/annotated/biotypes_go_raw.txt", index=None, sep="\\t")
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
            bio_go=bio_go.sort_values(by=["ensembl_gene_id"],ascending=True)
            bio_go.to_csv("{{project_folder}}/deseq2_output/annotated/biotypes_go.txt", sep= "\\t", index=None)
            bio_go.to_csv("{{biotypes_go}}", sep= "\\t", index=None)

        except Exception as e:
            print("An error occurred:", e,"\\ncontinuing without biomart go annotations.")
            bio_go=pd.DataFrame(columns=["ensembl_gene_id"])         

    else:
        shutil.copy( "{{biotypes_go}}".replace("biotypes_go.txt", "biotypes_go_raw_topgo.txt"), "{{project_folder}}/deseq2_output/annotated/biotypes_go_raw_topgo.txt" )
        shutil.copy( "{{biotypes_go}}".replace("biotypes_go.txt", "biotypes_go_raw.txt"), "{{project_folder}}/deseq2_output/annotated/biotypes_go_raw.txt" )
        shutil.copy( "{{biotypes_go}}", "{{project_folder}}/deseq2_output/annotated/biotypes_go.txt")
""",
    desc={
        "gtf":"",
        "biomart_dataset": "",
        "biomart_host": "",
        "project_folder":"",
        "biotypes_go":"'<path/to/biotypes_go.txt' inclusive 'biotypes_go.txt'"
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "1:00:00" }
)


parse_submission=jawm.Process(
    name="parse_submission",
    when=lambda p: not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "models.txt") ) ,
    script="""\
#!/usr/local/bin/python
import pandas as pd
from biomart import BiomartServer
import itertools
import AGEpy as age
import os
if not os.path.isdir("{{project_folder}}/deseq2_output/"):
    os.makedirs("{{project_folder}}/deseq2_output/")
sfile="{{samplestable}}"
if os.path.isfile(sfile):
    sdf=pd.read_excel(sfile)

elif os.path.isdir("{{kallisto_output}}"):
    print("Could not find a samples table of the form `pd.DataFrame(columns=['Files/Folders','group'])`")
    print("Workinf on predefined / expected kallisto output to generate test comparisons.")
    folders=os.listdir("{{kallisto_output}}")
    folders=[ s for s in folders if os.path.isdir(f"{{kallisto_output}}/{s}") ]
    folders=[ s for s in folders if not s.startswith("tmp.") ]
    folders=[ s for s in folders if ".Rep_" in s ]
    if not folders :
        raise Exception("Could not find a samples table of the form `pd.DataFrame(columns=['Files/Folders','group'])` nor a kallisto folder with Replicate information eg. <group>.Rep_<number> in it's sub folder names")
    folders.sort()
    groups=[ s.split(".Rep_")[0] for s in folders ]
    sdf=pd.DataFrame( { "Files/Folders":folders , "group":groups } )

sam_df=sdf.copy()
files_col=sam_df.columns.tolist()[0]
cond_col=sam_df.columns.tolist()[1]
sam_df.index=[s.split("{{read1_suffix}}")[0] for s in sam_df[files_col].tolist()]
sam_df=sam_df[[cond_col]]
sam_df.to_csv("{{project_folder}}/deseq2_output/samples_MasterTable.txt", sep="\t")
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
            print("{{project_folder}}/deseq2_output/"+filename+".input.tsv")
            outdf_.to_csv("{{project_folder}}/deseq2_output/"+filename+".input.tsv", sep="\t")
            text=[ filename, m, g, str(ref), coef ]
            text="\\t".join(text)
            textout.append(text)
        
with open("{{project_folder}}/deseq2_output/models.txt", "w") as mout:
    mout.write("\\n".join(textout) + "\\n")
""",
    desc={
        "samplestable":"",
        "project_folder":"",
        "read1_suffix":""
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "4GB", "-t": "1:00:00" }
)

tx2gene=jawm.Process(
    name="tx2gene",
    when=lambda p: not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "tx2gene.csv") ) ,
    script="""\
#!/usr/local/bin/python
import AGEpy as age
import os
if not os.path.isdir("{{project_folder}}/deseq2_output/"):
    os.makedirs("{{project_folder}}/deseq2_output/")
GTF=age.readGTF("{{gtf}}")
GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
GTF["transcript_id"]=age.retrieve_GTF_field(field="transcript_id",gtf=GTF)
tx2gene=GTF[["transcript_id","gene_id"]].drop_duplicates().dropna()
tx2gene.columns=["TXNAME","GENEID"]
tx2gene[["TXNAME","GENEID"]].to_csv("{{project_folder}}/deseq2_output/tx2gene.csv", quoting=1, index=None)
""",
    desc={
        "gtf":"",
        "project_folder":"", 
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "4GB", "-t": "1:00:00" }
)


tx2gene_proc=jawm.Process(
    name="tx2gene_proc",
    when=lambda p: not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "deseq2.part1.Rdata") ) ,
    script="""\
#!/usr/bin/Rscript
library(tidyverse)
library(tximportData)
library(tximport)
library(DESeq2)
library(rhdf5)
library(readr)
library(apeglm)
tx2gene <- read_csv("{{project_folder}}/deseq2_output/tx2gene.csv")
tx2gene <- tx2gene[rowSums(is.na(tx2gene)) == 0,]
tx2gene <- droplevels(tx2gene)
save.image("{{project_folder}}/deseq2_output/deseq2.part1.Rdata")
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
    when=lambda p: not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , p.var["input_file"].split(".input.tsv")[0]+".results.tsv" ) ) ,
    script="""\
#!/usr/bin/Rscript
load("{{project_folder}}/deseq2_output/deseq2.part1.Rdata")
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

out=paste("{{project_folder}}/deseq2_output/",coef,".results.tsv", sep="")

# filein
sampleTable<-read.delim2("{{project_folder}}/deseq2_output/{{input_file}}",sep = "\t", row.names = 1)
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
  circRNA = read.delim('{{project_folder}}/{{circRNA}}/CircRNACount', check.names = FALSE, as.is = TRUE )
  coord = read.delim('{{project_folder}}/{{circRNA}}/CircCoordinates', check.names = FALSE, as.is = TRUE)
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
    desc={
        "input_file":"",
        "project_folder":"",
        "kallisto_output":"",
        "circRNA": ""
    },
    container="mpgagebioinformatics/deseq2:1.38.0",
    manager_slurm={ "-c": 20, "--mem": "20GB", "-t": "1:00:00" }
)

mastertable=jawm.Process(
    name="mastertable",
    when=lambda p: not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "all_results_stats.xlsx" ) ) ,
    script="""\
#!/usr/bin/Rscript
load("{{project_folder}}/deseq2_output/deseq2.part1.Rdata")
library(tidyverse)
library(tximportData)
library(tximport)
library(DESeq2)
library(rhdf5)
library(readr)
library(apeglm)

sampleTable<-read.delim2("{{project_folder}}/deseq2_output/samples_MasterTable.txt",sep = "\\t", row.names = 1)
samples<-row.names(sampleTable)
# dir
dir <- "{{kallisto_output}}"
files<-file.path(dir, samples, "abundance.tsv")
names(files)<-samples
txi <- tximport(files, type = "kallisto", txOut = TRUE)
txi <- summarizeToGene(txi, tx2gene = tx2gene)
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
  circRNA = read.delim('{{project_folder}}/{{circRNA}}/CircRNACount', check.names = FALSE, as.is = TRUE )
  coord = read.delim('{{project_folder}}/{{circRNA}}/CircCoordinates', check.names = FALSE, as.is = TRUE)

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
openxlsx::write.xlsx(res_counts, "{{project_folder}}/deseq2_output/all_res_counts.xlsx", row.names = TRUE, col.names = TRUE)
write.table(res_counts, "{{project_folder}}/deseq2_output/all_res_counts.tsv", sep = "\\t", quote = T, row.names = T)
result_tables <- list.files('{{project_folder}}/deseq2_output/', pattern = '.results.tsv')
for(f in result_tables){
  tmp <- read.delim(paste0('{{project_folder}}/deseq2_output/', f))
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
res_counts = res_counts[,-1]
write.table(res_counts, "{{project_folder}}/deseq2_output/all_results_stats.tsv", sep = "\\t", quote = F, row.names = F)
openxlsx::write.xlsx(res_counts, "{{project_folder}}/deseq2_output/all_results_stats.xlsx", row.names = FALSE, col.names = TRUE)

sessionInfo()
""",
    desc={
        "circRNA": "",
        "project_folder":"",
        "kallisto_output":"",
    },
    container="mpgagebioinformatics/deseq2:1.38.0",
    manager_slurm={ "-c": 8, "--mem": "40GB", "-t": "4:00:00" }
)

annotator=jawm.Process(
    name="annotator",
    when=lambda p: not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "annotated", "masterTable_annotated.xlsx" ) ) ,
    script="""\
#!/usr/local/bin/python
import pandas as pd
import os
import AGEpy as age
import shutil

if not os.path.isdir("{{project_folder}}/deseq2_output/annotated/") :
    os.makedirs("{{project_folder}}/deseq2_output/annotated/")

if os.path.exists("{{project_folder}}/deseq2_output/annotated/biotypes_go.txt"):
    bio_go=pd.read_csv( "{{project_folder}}/deseq2_output/annotated/biotypes_go.txt", sep="\t")
else:
    bio_go=pd.DataFrame(columns=["ensembl_gene_id"])

GTF=age.readGTF("{{gtf}}")
GTF["gene_id"]=age.retrieve_GTF_field(field="gene_id",gtf=GTF)
GTF["gene_name"]=age.retrieve_GTF_field(field="gene_name",gtf=GTF)
id_name=GTF[["gene_id","gene_name"]].drop_duplicates()
id_name.reset_index(inplace=True, drop=True)
id_name.columns=["ensembl_gene_id","gene_name"]
#id_name=pd.read_table("{{project_folder}}/kallisto_index/cdna.norRNA.tsv")
#id_name=id_name[["gene_id","gene_symbol","description"]]
#id_name.columns=["ensembl_gene_id","gene_name","description"]
#id_name=id_name.drop_duplicates()
#id_name.reset_index(inplace=True,drop=True)

deg_files=os.listdir("{{project_folder}}/deseq2_output/")
deg_files=[ s for s in deg_files if "results.tsv" in s ]
i=1
s=[]
dfs={}
for f in deg_files:
    df=pd.read_table("{{project_folder}}/deseq2_output/"+f)
    df=pd.merge(id_name,df,left_on=["ensembl_gene_id"],right_index=True, how="right") # change to gene_id
    df=pd.merge(df,bio_go,on=["ensembl_gene_id"],how="left")
    df=df.sort_values(by=["padj"],ascending=True)
    df.to_csv("{{project_folder}}/deseq2_output/annotated/"+f, sep="\\t",index=None)
    df.to_excel("{{project_folder}}/deseq2_output/annotated/"+f.replace('.tsv', '.xlsx'), index=None)
    n=f.split(".results.tsv")[0]
    s.append([i,n])
    df=df[df["padj"]<0.05]
    df.reset_index(inplace=True, drop=True)
    dfs[i]=df
    i=i+1
sdf=pd.DataFrame(s,columns=["sheet","comparison"])
EXC=pd.ExcelWriter("{{project_folder}}/deseq2_output/annotated/significant.xlsx")
sdf.to_excel(EXC,"summary",index=None)
for k in list(dfs.keys()):
    dfs[k].to_excel(EXC, str(k),index=None)
EXC.close()
mt=pd.read_csv("{{project_folder}}/deseq2_output/all_results_stats.tsv", sep="\\t")
mt_ann=pd.merge(id_name,mt,on=["ensembl_gene_id"], how="right")
mt_ann=pd.merge(mt_ann,bio_go,on=["ensembl_gene_id"],how="left")
mt_ann.to_csv("{{project_folder}}/deseq2_output/annotated/masterTable_annotated.tsv", sep="\\t",index=None)
mt_ann.to_excel("{{project_folder}}/deseq2_output/annotated/masterTable_annotated.xlsx", index=None)
""",
    desc={
        "gtf":"",
        "biomart_host": "",
        "project_folder":"",
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 2, "--mem": "4GB", "-t": "1:00:00" }
)

david=jawm.Process(
    name="david",
    when=lambda p: (  ( not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "annotated", "david.touch" ) ) )  &  ( p.var["DAVIDUSER"] != "" ) ) ,
    script="""\
#!/usr/local/bin/python
import pandas as pd
import AGEpy as age
import os 
import sys
from pathlib import Path
file_path = Path("{{project_folder}}/deseq2_output/annotated/david.touch")

deseq2="{{project_folder}}/deseq2_output/annotated/"
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
    desc={
        "DAVIDUSER": "",
        "daviddatabase": "",
        "project_folder": ""
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "12:00:00" }
)

topgo=jawm.Process(
    name="topgo",
    when=lambda p: not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "annotated", p.var["input_file"].replace("results.tsv","topGO.tsv") ) ) ,
    script="""\
#!/usr/bin/Rscript
library(topGO)
library(biomaRt)
library(plyr)
library(openxlsx)
## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fin="{{project_folder}}/deseq2_output/annotated/{{input_file}}" # topGO.tsv
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
    goterms<-read.delim("{{project_folder}}/deseq2_output/annotated/biotypes_go_raw_topgo.txt", as.is = TRUE)
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
        D_out[,'percent'] = NA
        D_out[, 'listTotals'] = length(sigGenes)
        D_out[, 'listTotals_used'] = length(sigGenes(sampleGOdata))
        D_out[, 'popTotals'] = length(allGenes(sampleGOdata))
        D_out[, 'popTotals_used'] = numGenes(sampleGOdata)
        D_out[, 'foldEnrichment'] = (D_out[,'Significant']/D_out[,'listTotals_used']) / (D_out[, 'Annotated']/D_out[, 'popTotals_used'])
        D_out[, 'bonferroni'] = p.adjust(D_out[,'classicFisher'], method = "bonferroni")
        D_out[, 'benjamini'] = p.adjust(D_out[,'classicFisher'], method = "BY")
        D_out[, 'afdr'] = p.adjust(D_out[,'classicFisher'], method = "fdr")


        ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        D_out = D_out[, c('categoryName', 'termName', 'Significant', 'percent', 'classicFisher', 'geneIds',
                            'listTotals', 'listTotals_used', 'Annotated', 'popTotals', 'popTotals_used', 
                            'Expected',  'foldEnrichment', 'bonferroni', 'benjamini', 'afdr', "gn", "gfc")]

        names(D_out) <- c('categoryName', 'termName', 'listHits', 'percent', 'classicFisher', 'geneIds',
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
        "project_folder":"",
        "input_file":"",
    },
    container="mpgagebioinformatics/topgo:2.50.0",
    manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "4:00:00" }
)

cellplot=jawm.Process(
    name="cellplot",
    when=lambda p: ( not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "annotated", p.var["inFile"].split( p.var["filetype"] )[0] )+p.var["category"]+".cellplot.pdf"  ) ) & \
                   ( not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "annotated", p.var["inFile"].split( p.var["filetype"] )[0] )+p.var["category"]+".cellplot.touch" ) ) ,
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

setwd("{{project_folder}}/deseq2_output/annotated/")

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
  D$ease <- as.numeric(as.character(D$classicFisher))
}

D$ease <- as.numeric(as.character(D$ease))
D$foldEnrichment <- as.numeric(as.character(D$foldEnrichment))
D$listHits <- as.numeric(as.character(D$listHits))

# Added for handling cellplot NA values
# Handle log2fc properly
D$log2fc <- gsub("inf", "Inf", as.character(D$log2fc))
D$log2fc <- as.numeric(D$log2fc)

# Remove rows where `ease`, `foldEnrichment`, or `log2fc` are NA
D <- D[!is.na(D$ease) & !is.na(D$foldEnrichment) & !is.na(D$log2fc), ]

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
        "project_folder":"",
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
    when=lambda p: ( not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "annotated", p.var["input_file"].replace( ".results.tsv", ".RcisTarget.xlsx")  ) ) &  ( p.var["rcis_db"] != "" ) ),
    script="""\
#!/usr/bin/Rscript
library(RcisTarget)
library(openxlsx)
setwd("{{project_folder}}/deseq2_output/annotated/")
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
    desc={
        "rcis_db": "Available from https://resources.aertslab.org/cistarget/",
        "project_folder":"",
        "input_file":""
    },
    container="mpgagebioinformatics/rcistarget:1.17.0",
    manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "4:00:00" }
)

qc_plots=jawm.Process(
    name="qc_plots",
    when=lambda p: not os.path.isfile( os.path.join( p.var["project_folder"], "qc_plots" , "pca_all_samples.pdf"  ) ),
    script="""\
#!/bin/bash
mkdir -p {{project_folder}}/qc_plots/
QC_plots -o {{project_folder}}/qc_plots/ -de {{project_folder}}/deseq2_output/ -s {{project_folder}}/deseq2_output/samples_MasterTable.txt -t {{project_folder}} -sp {{spec}}
""",
    desc={
        "project_folder":"",
        "spec": ""
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "8GB", "-t": "4:00:00" }
)

get_ip=jawm.Process(
    name="get_ip",
    when=lambda p: ( not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "annotated", "string.done"  ) ) &  ( p.var["cytoscape_host"] != "" ) ) ,
    script="""\
#!/bin/bash
while [[ ! -f {{cytoscape_host}} ]] ; do 
  echo "waiting for cytoscape to be available"
  sleep 3$((RANDOM % 9))
done
mv {{cytoscape_host}} {{cytoscape_host}}_inuse
touch {{project_folder}}/deseq2_output/annotated/string.running
""",
    desc={
        "cytoscape_host":""
    },
)

# using string over cytsocape is deprecated 
# a new module without cytoscape is currently under development

string=jawm.Process(
    name="string",
    when=lambda p: ( not os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "annotated", "string.done"  ) ) &  ( p.var["cytoscape_host"] != "" ) ),
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
    input_files=os.listdir("{{project_folder}}/deseq2_output/annotated/")
    input_files=[s for s in input_files if ".results.tsv" in s ]
    input_files=[ os.path.join("{{project_folder}}/deseq2_output/annotated/",s) for s in input_files if ".results.tsv" in s ]
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
    input_files=os.listdir("{{project_folder}}/deseq2_output/annotated/")
    input_files=[s for s in input_files if ".results.tsv" in s ]
    input_files=[ os.path.join("{{project_folder}}/deseq2_output/annotated/",s) for s in input_files if ".results.tsv" in s ]
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

file_path = Path("{{project_folder}}/deseq2_output/annotated/string.done" )
file_path.touch()
""",
    desc={
        "cytoscape_host":"Path to a txt file containing one line with the IP address of your running cytoscape instance.",
        "biomart_host": "",
        "species": "",
        "project_folder":""
    },
    container="mpgagebioinformatics/rnaseq.python:3.8-8",
    manager_slurm={ "-c": 1, "--mem": "20GB", "-t": "4:00:00" }
)

release_ip=jawm.Process(
    name="release_ip",
    when=lambda p: ( os.path.isfile( os.path.join( p.var["project_folder"], "deseq2_output" , "annotated", "string.running"  ) ) &  ( p.var["cytoscape_host"] != "" ) ),
    script="""\
#!/bin/bash
if [[ -f {{cytoscape_host}}_inuse ]] ; then mv {{cytoscape_host}}_inuse {{cytoscape_host}} ; fi
rm -rf {{project_folder}}/deseq2_output/annotated/string.running
""",
    desc={
        "cytoscape_host": "",
        "project_folder":""
    },
    container=""
)

upload_paths=jawm.Process(
    name="upload_paths",
    when=True,
    script="""\
rm -rf upload.txt

cd {{project_folder}}/deseq2_output/annotated

for f in $(ls *.results.xlsx) ; do echo "deseq2 $(readlink -f ${f})" >>  upload.txt_ ; done
echo "deseq2 $(readlink -f significant.xlsx)" >>  upload.txt_
echo "deseq2 $(readlink -f masterTable_annotated.xlsx)" >>  upload.txt_

if [[ $(ls  | grep cytoscape) ]] ; then
  for f in $(ls *.cytoscape.* ) ; do echo "cytoscape $(readlink -f ${f})" >>  upload.txt_ ; done
fi

if [[ $(ls  | grep DAVID) ]] ; then
  for f in $(ls *.DAVID.* ) ; do echo "david $(readlink -f ${f})" >>  upload.txt_ ; done
fi

if [[ $(ls  | grep RcisTarget) ]] ; then
  for f in $(ls *.RcisTarget.* ) ; do echo "rcistarget $(readlink -f ${f})" >>  upload.txt_ ; done
fi

if [[ $(ls  | grep topGO) ]] ; then
  for f in $(ls *.topGO.* ) ; do echo "topgo $(readlink -f ${f})" >>  upload.txt_ ; done
fi

uniq upload.txt_ upload.txt 
rm upload.txt_

cd {{project_folder}}/qc_plots
rm -rf upload.txt 
for f in $(ls *.* | grep -v upload.txt) ; do echo "qc_plots $(readlink -f ${f})" >>  upload.txt_ ; done

uniq upload.txt_ upload.txt 
rm upload.txt_
""",
    desc={
        "project_folder": ""
    },
    container=""
)

if __name__ == "__main__":
    import sys
    from jawm.utils import workflow
    from pathlib import Path
    import glob

    workflows, var, args, unknown_args = jawm.utils.parse_arguments(['main','deseq2','test', "cistarget"])

    if workflow(["main","deseq2","test","cistarget"], workflows):

        tx2gene.execute( )
        annotations.execute()
        parse_submission.execute()
        tx2gene_proc.execute(tx2gene.hash)

        jawm.Process.wait([ tx2gene_proc.hash, parse_submission.hash ])

        base = os.path.join(deseq2.var["project_folder"], "deseq2_output")
        tests = glob.glob(os.path.join(base, "*input.tsv"))

        deseq2_jobs=[]

        for file in tests:

            # clone the processes required for each file
            deseq2_=deseq2.clone()

            deseq2_.var["input_file"]=os.path.basename(file)
            deseq2_.execute()

            deseq2_jobs.append( deseq2_.hash )
        
        mastertable.execute( )
        deseq2_jobs.append( mastertable.hash )

        annotator.execute( deseq2_jobs )


        # wait for annotator to complete before listing input files and starting david
        # as they are generated by annotator
        jawm.Process.wait( annotator.hash )

        david.execute( )

        base = os.path.join(deseq2.var["project_folder"], "deseq2_output", "annotated")
        tests_results = glob.glob(os.path.join(base, "*results.tsv"))
        topgo_jobs=[]
        for file in tests_results:
            topgo_=topgo.clone()
            topgo_.var["input_file"]=os.path.basename(file)
            topgo_.execute( )
            topgo_jobs.append( topgo_.hash )

        jawm.Process.wait( topgo_jobs + [ david.hash ] )

        base = os.path.join(deseq2.var["project_folder"], "deseq2_output", "annotated")
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

    if workflow("test", workflows):

        with open(os.path.join(var["project_folder"], "test.txt"), 'w') as out:
            out.write("Test completed.")

        # for the test workflow we might also do something more
        print("Test completed.")
