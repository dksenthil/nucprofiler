## NUCPROFILER:

DNA shape profile for genomic fragments.

NUCPROFILER has two scripts. 

- **nucprofiler_intraparams.R** calcuates intra base pair parameters like  "ProT","Stretch","Buckle","Shear","Opening","Stagger" as well as Electrosatic potnetial (EP) and minor groove width (MGW).

- **nucprofiler_interbaseparams.R** calcuates intra base pair parameters like "HelT","Roll","Tilt","Shiftr","Slide", and "Rise".




## Prerequisite R packages
```
library("seqinr")
library(stringi)
library(tools)
library(data.table)

```
## Usage

```

user@machine:~nucprofiler$ Rscript nucprofiler_intraparams.R      example/small.fasta 
user@machine:~nucprofiler$ Rscript nucprofiler_interbaseparams.R  example/small.fasta

```
## A sample Bash script (batch)

An bash script to run in batch is provied.

```
#!/bin/bash
for file in `cat list_of_fasta`
do
  
  Rscript nucprofiler_interbaseparams.R $file
  Rscript nucprofiler_intraparams.R     $file
  
done

```

## Query table information 

The query table is inside *data* folder.


## References

*  J. Li et al. Expanding the repertoire of DNA shape features for  genome-scale studies of transcription factor binding. Nucleic Acids Res. 45, 12877-12887 (2017)

* T.P. Chiu et al. Genome-wide prediction of minor-groove electrostatic potential enables biophysical modeling of protein-DNA binding.Nucleic Acids Res. 45, 12565-12576 (2017)

## Example

The *example* folder has there types of files representative of small, medium and big fasta files.

## Sample R Markdown documents

For nucprofiler_intraparams.R the R makrdown is provided.

## New in the latest version

Version: 1.0.1
Date: 2018-06-01

* fixed some typos and paths
* added query table publication
* added example fasta files

## Outlook

1. Add Matrix plotting option.
2. Add dinucleotide based porperty.