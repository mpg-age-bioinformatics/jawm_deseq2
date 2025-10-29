# jawm_deseq2

This is a jawm deseq2 module.

Installing jawm:
```
pip install git+ssh://git@github.com/mpg-age-bioinformatics/jawm.git
```
For more information on jawm please visit jawm's repo on [GitHub.com](https://github.com/mpg-age-bioinformatics/jawm/tree/main).

Example usage:
```
# clone this module
git clone git@github.com:mpg-age-bioinformatics/jawm_deseq2.git

cd jawm_deseq2

# download test data
jawm-test -r download

# docker
jawm deseq2.py deseq2 -p ./yaml/docker.yaml

# slurm & apptainer with multiple yaml files
jawm deseq2.py deseq2 -p ./yaml/vars.yaml ./yaml/hpc.yaml
```

Testing this module on your system's python, jawm, and docker installations:
```
jawm-test --python_versions system --jawm_versions local
```
More information on running and developing tests can be found in `./test/README.md`.

Additional jawm workflows are available [here (GitHub.com)](https://github.com/mpg-age-bioinformatics?q=jawm_&type=all&language=&sort=).
