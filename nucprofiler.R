## Loop throught miltifasta
library("seqinr")
library(stringi)
library(tools)
library(data.table)



## Get the input file as commandline argument and read the file
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop ("Usage: Rscipt  run_DNAshapeProfiler.R <Fastafile> ")
inputfile <- args[1]

#inputfile <-c("temp.fasta")

fastafile <- read.fasta(
  file = inputfile,
  as.string = TRUE,
  set.attributes = TRUE,
  strip.desc = FALSE,
  seqonly = FALSE
)


## Read the DNAshape query table
rohdata <-
  fread("../data/DNA_shape_query_table.csv")
framedata <- data.frame(rohdata)
rownames(framedata) <- framedata[, 1]
proptable <- framedata[1:512, ]


## Define function for handling reverse complementary pentamer
revcomp <- function(nucSeq)
  return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq)))


## DNAshape_retreive function1:
dnashape_shapemat1 <- function(pentamer_vec, propcol) {
  shape_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  part2 <- c(NA)
  dim(part2) <- c(1, 1)
  for (b in 1:length(pentamer_vec)) {
    checkvar <- grepl(pentamer_vec[b], rownames(proptable))
    
    if (any(checkvar)) {
      part1 <-
        proptable[grepl(pentamer_vec[b], rownames(proptable)),][propcol]
      temp <- c(part1[1, 1], part2[1, 1])
      shape_mat[b] <- mean(temp, na.rm = TRUE)
      part2 <-
        proptable[grepl(pentamer_vec[b], rownames(proptable)),][propcol +
                                                                  1]
      if(b == length(pentamer_vec)){shape_mat[b+1] <- part2}
    } else {
      part1 <-
        proptable[grepl(revcomp(pentamer_vec[b]), rownames(proptable)),][propcol +
                                                                           1]
      temp <- c(part1[1, 1], part2[1, 1])
      shape_mat[b] <- mean(temp, na.rm = TRUE)
      part2 <-
        proptable[grepl(revcomp(pentamer_vec[b]), rownames(proptable)),][propcol]
      
      if (b == length(pentamer_vec)) {
        shape_mat[b + 1] <- part2
      }
    }
  }
  return(shape_mat)
}
###--------
## DNAshape_retreive function2:
dnashape_shapemat2 <- function(pentamer_vec, propcol) {
  shape_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  for (b in 1:length(pentamer_vec)) {
    checkvar <- grepl(pentamer_vec[b], rownames(proptable))
    
    if (any(checkvar)) {
      shape_mat[b] <-
        proptable[grepl(pentamer_vec[b], rownames(proptable)),][propcol]
    } else {
      shape_mat[b] <-
        proptable[grepl(revcomp(pentamer_vec[b]), rownames(proptable)),][propcol]
                                                            
    }
  }
  
  return(c(NA, shape_mat, NA))
}

## Main script

 rownum <- c(4,6,9,14,16,19,2,3,8,11,12,13,18,21)
 namenum <- c("Roll","HelT","Tilt","Rise","Shift","Slide","MGW","ProT","Stretch","Buckle","Shear","Opening","Stagger","EP")

for(h in 1:length(rownum)){
  # remove files if exist in the cwd
  outname <- file_path_sans_ext(basename(inputfile))
  shapefilename <- paste(outname, ".", namenum[h], ".mx", sep = "")
  if (file.exists(shapefilename)){
    file.remove(shapefilename)
  }
array1 <- vector()
for (a in 1:length(fastafile)) {
  print(names(fastafile)[a])
  pentamer_vec <- vector()
  sequence <- paste(toupper(fastafile[a]))
  seqvec <- strsplit(sequence, "")
  
  for (i in 1:(nchar(sequence) - 4)) {
    pentamer_vec[i] <- paste0(seqvec[[1]][i:(i + 4)], collapse = "")
  }
  headname <-
    array(c(paste0(getAnnot(fastafile)[a])), dim = c(1, 1))
  #shapeout <- dnashape_shapemat1(pentamer_vec, rownum[h])
  if ((rownum[h] == 4 ) || (rownum[h] == 6 ) || (rownum[h] == 9 ) || (rownum[h] == 14 ) ||(rownum[h] == 16 ) || (rownum[h] == 19 )){
    shapeout <- dnashape_shapemat1(pentamer_vec, rownum[h])
  }else if ((rownum[h] == 2 ) || (rownum[h] == 3 ) || (rownum[h] == 8 ) || (rownum[h] == 11 ) ||(rownum[h] == 12 ) || (rownum[h] == 13 ) || (rownum[h] == 18 ) || (rownum[h] == 21 ) ){
    shapeout <- dnashape_shapemat2(pentamer_vec, rownum[h])
  }else {print ("Error in input!","\n")}
  
  array1 <- paste0(c(headname[1], "NA", round(as.numeric(rbind(shapeout)), digits =2 ), "NA"), sep = "")
  write.table( rbind(array1), file = shapefilename, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE, na = "NA",  dec = ".", quote = FALSE )
  
}# fasta loop close

}

