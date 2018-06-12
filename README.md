# nucprofiler

# MDplot: Visualise Molecular Dynamics Analyses
NUCPROFILER has two set of script. 

nucprofiler_intraparams.R calcuates intra base pair parameters like 
"ProT","Stretch","Buckle","Shear","Opening","Stagger" as well as Electrosatic potnetial (EP) and minor groove width (MGW).



```
user@machine:~nucprofiler$ Rscript nucprofiler_intraparams.R example/temp.fasta 

user@machine$ Rscript nucprofiler_intrabaseparams.R filename.fasta

```

## Call from within bash script
An Rscript interface is provided, allowing to set most options:
```
#!/bin/bash
for file in `cat list_of_fasta`
do
  
  Rscript nucprofiler_intrabaseparams.R $file
  Rscript nucprofiler_intraparams.R     $file
done

```

## Additional information and examples
The query table is in data folder

The example folder has there types of files representative of samll. medium and big fasta files.


## New in the latest major version
Version: 1.0.1
Date: 2017-07-04

* fixed some typos in the manual pages
* added vignette (publication)
* added special input support for function "load_timeseries()" to be able
  to load multi-column timeseries data

Version: 1.0.0
Date: 2017-02-24

* fixed issue with proper residue display when sub-selection was done in function 'dssp()'
* fixed issue with bin-expansion in function 'load_noe()'
* removed unnecessary input parameters from several functions
* added 'stride' to function 'load_dssp_ts()'
* changed parameters of function 'load_rmsf()'