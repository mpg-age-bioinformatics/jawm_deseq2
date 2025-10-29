import jawm

deseq2_p1=jawm.Process( 
    name="deseq2_p1",
    script="""#!/bin/bash
echo "{{extra_args}} {{my_deseq2_argument}}" 2>&1 | tee {{output}}/deseq2.txt
""",

    # arguments for the script above :
    
    # var={
    #     "extra_args": "",
    #     "my_deseq2_argument":"This is just a deseq2.", 
    #     "mk.output":"<output_folder>", # the prefix "mk." leads to the creation of this folder and volume mapping if you are using containers
    # },

    # here you can describe your variables
    desc={
        "extra_args": "use this if you want to add not preset arguments",
        "my_deseq2_argument":"Some text that will be printed to the screen.", 
        "output":"Folder for output", # the prefix "mk." leads to the creation of this folder and volume mapping if you are using containers
    },
    
    # example arguments for jawn

    # manager="slurm",
    # manager_slurm={
    #     "-p":"cluster,dedicated", 
    #     "--mem":"20GB", 
    #     "-t":"1:00:00", 
    #     "-c":"8" 
    # },
    
    # container="docker://mpgagebioinformatics/fastqc:0.11.9",
    # environmnent="apptainer",
    # environment_apptainer={ '-B': [input_file, output_folder] }
    
    # container="mpgagebioinformatics/fastqc:0.11.9",
    # environmnent="docker",
    # environment_docker={ '-v': [input_file, output_folder] },


    # param_file="yaml/apptainer.params.yaml" ,
    # param_file=[ "yaml/apptainer.params.yaml" , "yaml/slurm.params.yaml" ],
  
)

deseq2_p2=jawm.Process( 
    name="deseq2_p2",
    script="""#!/usr/bin/env python3
with open("{{map.file}}", "r") as src, open("{{output}}/deseq2.txt", "a") as dst:
    dst.write(src.read())
""",

    # arguments for the script above :
    
    # var={
    #     "mk.output":"<output_folder>", # the prefix "mk." leads to the creation of this folder and volume mapping if you are using containers
    #     "map.file":"<some_file>"# the prefix "map." leads to the mapping of this file if you are using containers
    # },
  
)

deseq2_p3=jawm.Process( 
    name="deseq2_p3",
    script="""#!/usr/bin/env Rscript
write( "\nDemo completed", file = "{{output}}/deseq2.txt", append = TRUE)
""",

    # arguments for the script above :
    
    # var={
    #     "mk.output":"<output_folder>", # the prefix "mk." leads to the creation of this folder and volume mapping if you are using containers
    # },
  
)


if __name__ == "__main__":
    import sys
    from jawm.utils import workflow
    from jawm.utils import load_modules

    # load modules from local folders
    load_modules("submodules")

    # load modules from online git repos
    load_modules("jawm_template")

    # it can be used in the form
    # load_modules(["modules","jawm_template@<tag/commit_hash>"])
    
    # or for latest available tag
    # load_modules("jawm_template@latest")

    workflows, var, args, unknown_args= jawm.utils.parse_arguments(["main","deseq2","test"],)

    # usage: 

    if workflow( ["main","deseq2","test"], workflows ) :

        # execute process
        deseq2_p1.execute()

        # execute a process with dependencies
        deseq2_p2.depends_on=[deseq2_p1.hash]
        deseq2_p2.execute()

        # wait for all above processes to complete
        jawm.Process.wait(deseq2_p2.hash)

        # print the output
        print(deseq2_p1.get_output())
        print(deseq2_p2.get_output())

        deseq2_submodule.deseq2_submodule_p1.execute()
        template._template_p2.execute()

        jawm.Process.wait([ deseq2_submodule.deseq2_submodule_p1.hash, template._template_p2.hash  ])
        print(deseq2_submodule.deseq2_submodule_p1.get_output())
        print(template._template_p2.get_output())

    if workflow( "test", workflows ) :

        # for the test workflow we also do something more (just for deseq2)
        deseq2_p3.execute()
        jawm.Process.wait( deseq2_p3.hash)
        print("Test completed.")


    sys.exit(0)
