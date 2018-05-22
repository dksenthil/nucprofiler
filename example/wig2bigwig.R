## convert WIG to BIGWIG
library(rtracklayer)
library(GenomeInfoDb) 
# define human genome source
HG19 <- Seqinfo(genome="hg19")

Wiglist <- list.files(path = ".", pattern = ".wig", all.files = FALSE,
                                full.names = FALSE, recursive = FALSE,
                                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
# the folder where the results will be stored
out_folder=paste0("./out")
dir.create( out_folder, showWarnings = T, recursive = F )

for( h in 1:length(Wiglist)){
  newname <- file_path_sans_ext(Wiglist[h])
  NN <- paste0("./out/",newname,".bw")
  wigToBigWig(Wiglist[h], seqinfo = HG19, dest = NN, clip = TRUE)
}
