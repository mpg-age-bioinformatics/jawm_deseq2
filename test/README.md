# jawm tests

It is good working practices to prepare data driven tests of your workflows. 
This will allow you to continuously test your code against the lalest jawm 
version and at the same time give you a framework to test how changes to your 
code or dependent docker images might affect your results.

The file `tests.txt` lists all tests that should be performed:
```
#jawm_file.py;workflow;parameters.file1.yaml,parameters.file2.yaml;"Test name";test_hash
deseq2.py;test;./test/yaml/test.yaml;"Main workflow test";ce67151f2a16c1ee96e65c5eec6e7b07288f003fb3affce6f08568cf11a1525f
```
The `test_hash` is created based on the set of files or folders you define as the golden reference for your worklow. 
These are defined in `yaml` eg. `./yaml/test.yaml` file under the `- scope: hash` :
```
- scope: hash
  include: ./test/test-output/deseq2.txt
  overwrite: true # overwrite any existing hash already present on the logs folder
```
Initally, on your first test, you won't have an hash, leave it blanck and let jwam fill the value for you.
```
jawm-test
``` 
If you need to overwrite the hash on teh `tests.txt` file you can:
```
jawm-test -o
```
Additionaly, these hashes are stored under `./logs/jwam_hashes/<module_name>_user_defined.history` and 
`./logs/jwam_hashes/<module_name>.hash` for the latest hash.

If you are running multiple variations of the same module you might want to change the path of the logs directory with 
```
jawm -l <path_to_new_logs_directory>
```

The file `data.txt` lists all the files required for your tests in the form:
```
<md5sums>  <file_name>  <url_for_file_download> 
```
These will be downloaded when running `jawm-test`.

You can also trigger a workflow test from using GitHub's client eg.:
```
gh workflow run "Test Module"
```
then check the run status with
```
gh run list --workflow="test.yaml"
```
and watch your remote test with
```
gh run watch <run-id>
```
