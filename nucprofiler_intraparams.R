# Nucprofiler to calculate intra base pair parameters for Genomic fragments

library("seqinr")
library(stringi)
library(tools)
library(data.table)

## Get the input file as commandline argument and read the file
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop ("Usage: Rscipt  nucprofiler_intraparams.R <Fastafile> ")
inputfile <- args[1]

#inputfile <- c("example/temp.fasta")

fastafile <- read.fasta(
  file = inputfile,
  as.string = TRUE,
  set.attributes = TRUE,
  strip.desc = FALSE,
  seqonly = FALSE
)


## Read the DNAshape query table
rohdata <-
  fread("data/DNA_shape_query_table.csv")
framedata <- data.frame(rohdata)
rownames(framedata) <- framedata[, 1]
proptable <- framedata[1:512, ]


## Define function for handling reverse complementary pentamer
revcomp <- function(nucSeq)
  return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq)))


#rownum <- c(4, 6, 9, 14, 16, 19, 2, 3, 8, 11, 12, 13, 18, 21)
## DNAshape_retreive function2:
dnashape_shapemat <- function(pentamer_vec) {
  ## Mapping HASh-KEY 2,3,8,11,12,13,18,21
  
  map <- new.env(hash = T, parent = emptyenv())
  pent5 <- rownames(framedata)
  petnew <- head(pent5, -4)
  num  <- 1:512
  for (i in seq_along(petnew))
  {
    print
    key <- petnew[i]
    map[[key]] <- num[i]
    # print(ls(map))
  }
  #"MGW","ProT","Stretch","Buckle","Shear","Opening","Stagger","EP"
  MGW_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  
  ProT_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  
  Stretch_mat <-
    array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  
  Buckle_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  
  Shear_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  
  Opening_mat <-
    array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  
  Stagger_mat <-
    array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  
  EP_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  
  
  for (b in seq_along(pentamer_vec)) {
    checkvar <- grepl(pentamer_vec[b], rownames(proptable))
    
    if (any(checkvar)) {
      MGW_mat[b]     <- proptable[map[[pentamer_vec[b]]],][2]
      ProT_mat[b]    <- proptable[map[[pentamer_vec[b]]],][3]
      Stretch_mat[b] <- proptable[map[[pentamer_vec[b]]],][8]
      Buckle_mat[b]  <- proptable[map[[pentamer_vec[b]]],][11]
      Shear_mat[b]   <- proptable[map[[pentamer_vec[b]]],][12]
      Opening_mat[b] <- proptable[map[[pentamer_vec[b]]],][13]
      Stagger_mat[b] <- proptable[map[[pentamer_vec[b]]],][18]
      EP_mat[b]      <- proptable[map[[pentamer_vec[b]]],][21]
      
    } else {
      MGW_mat[b] <- proptable[map[[revcomp(pentamer_vec[b])]],][2]
      ProT_mat[b] <- proptable[map[[revcomp(pentamer_vec[b])]],][3]
      Stretch_mat[b] <-
        proptable[map[[revcomp(pentamer_vec[b])]],][8]
      Buckle_mat[b] <-
        proptable[map[[revcomp(pentamer_vec[b])]],][11]
      Shear_mat[b] <-
        proptable[map[[revcomp(pentamer_vec[b])]],][12]
      Opening_mat[b] <-
        proptable[map[[revcomp(pentamer_vec[b])]],][13]
      Stagger_mat[b] <-
        proptable[map[[revcomp(pentamer_vec[b])]],][18]
      EP_mat[b] <- proptable[map[[revcomp(pentamer_vec[b])]],][21]
      
    }
  }
  
  newList <-
    list(
      "MGW" = c(NA, MGW_mat, NA),
      "ProT" = c(NA, ProT_mat, NA),
      "Stretch" = c(NA, Stretch_mat, NA),
      "Buckle" = c(NA, Buckle_mat, NA),
      "Shear" = c(NA, Shear_mat, NA),
      "Opening" = c(NA, Opening_mat, NA),
      "Stagger" = c(NA, Stagger_mat, NA),
      "EP" = c(NA, EP_mat, NA)
    )
  
  
  return(newList)
}

## Main script
# remove files if exist in the cwd
outname <- file_path_sans_ext(basename(inputfile))

mgwfilename <-
  paste(outname, ".", "MGW.mx", sep = "")
mgwfilenameold <-
  paste(outname, ".", "MGW.mx.old", sep = "")
if (file.exists(mgwfilename)) {
  file.rename(mgwfilename, mgwfilenameold)
}
protfilename <-
  paste(outname, ".", "ProT.mx", sep = "")
protfilenameold <-
  paste(outname, ".", "ProT.mx.old", sep = "")
if (file.exists(protfilename)) {
  file.rename(protfilename, protfilenameold)
}
stretchfilename <-
  paste(outname, ".", "Stretch.mx", sep = "")
stretchfilenameold <-
  paste(outname, ".", "Stretch.mx.old", sep = "")
if (file.exists(stretchfilename)) {
  file.rename(stretchfilename, stretchfilenameold)
}
bucklefilename <-
  paste(outname, ".", "Buckle.mx", sep = "")
bucklefilenameold <-
  paste(outname, ".", "Buckle.mx.old", sep = "")
if (file.exists(bucklefilename)) {
  file.rename(bucklefilename, bucklefilenameold)
}
shearfilename <-
  paste(outname, ".", "Shear.mx", sep = "")
shearfilenameold <-
  paste(outname, ".", "Shear.mx.old", sep = "")
if (file.exists(shearfilename)) {
  file.rename(shearfilename, shearfilenameold)
}
openingfilename <-
  paste(outname, ".", "Opening.mx", sep = "")
openingfilenameold <-
  paste(outname, ".", "Opening.mx.old", sep = "")
if (file.exists(openingfilename)) {
  file.rename(openingfilename, openingfilenameold)
}
staggerfilename <-
  paste(outname, ".", "Stagger.mx", sep = "")
staggerfilenameold <-
  paste(outname, ".", "Stagger.mx.old", sep = "")
if (file.exists(staggerfilename)) {
  file.rename(staggerfilename, staggerfilenameold)
}
epfilename <-
  paste(outname, ".", "EP.mx", sep = "")
epfilenameold <-
  paste(outname, ".", "EP.mx.old", sep = "")
if (file.exists(epfilename)) {
  file.rename(epfilename, epfilenameold)
}

for (a in seq_along(fastafile)) {
  print(names(fastafile)[a])
  pentamer_vec <- vector()
  sequence <- paste(toupper(fastafile[a]))
  seqvec <- strsplit(sequence, "")
  
  pentamer_vec = sapply(1:(nchar(sequence) - 4), function(i) {
    paste0(seqvec[[1]][i:(i + 4)], collapse = "")
  })
  
  
  headname <-
    array(c(paste0(getAnnot(fastafile)[a])), dim = c(1, 1))
  listshape <- dnashape_shapemat(pentamer_vec)
  
  mgwarr <-
    paste0(c(headname[1], "NA", round(as.numeric(listshape$MGW), digits = 2), "NA"), sep = "")
  write.table(
    rbind(mgwarr),
    file = mgwfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  protarr <-
    paste0(c(headname[1], "NA", round(as.numeric(listshape$ProT), digits = 2), "NA"), sep = "")
  write.table(
    rbind(protarr),
    file = protfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  stretcharr <-
    paste0(c(headname[1], "NA", round(
      as.numeric(listshape$Stretch), digits = 2
    ), "NA"), sep = "")
  write.table(
    rbind(stretcharr),
    file = stretchfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  bucklearr <-
    paste0(c(headname[1], "NA", round(as.numeric(listshape$Buckle), digits = 2), "NA"), sep = "")
  write.table(
    rbind(bucklearr),
    file = bucklefilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  sheararr <-
    paste0(c(headname[1], "NA", round(as.numeric(listshape$Shear), digits = 2), "NA"), sep = "")
  write.table(
    rbind(sheararr),
    file = shearfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  openingarr <-
    paste0(c(headname[1], "NA", round(
      as.numeric(listshape$Opening), digits = 2
    ), "NA"), sep = "")
  write.table(
    rbind(openingarr),
    file = openingfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  staggerarr <-
    paste0(c(headname[1], "NA", round(
      as.numeric(listshape$Stagger), digits = 2
    ), "NA"), sep = "")
  write.table(
    rbind(staggerarr),
    file = staggerfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  eparr <-
    paste0(c(headname[1], "NA", round(as.numeric(listshape$EP), digits = 2), "NA"), sep = "")
  write.table(
    rbind(eparr),
    file = epfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
}

