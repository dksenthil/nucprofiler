#!/bin/bash
for file in `cat list_of_fasta`
do
  
  Rscript nucprofiler_interbaseparams.R $file
  Rscript nucprofiler_intraparams.R     $file
done