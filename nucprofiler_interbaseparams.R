## Nucprofiler to calculate inter base pair parameters for Genomic fragments

library("seqinr")
library(stringi)
library(tools)
library(data.table)

# # Get the input file as commandline argument and read the file
#  args <- commandArgs(trailingOnly = TRUE)
#  if (length(args) < 1)
#    stop ("Usage: Rscipt  nucprofiler.R <Fastafile> ")
#  inputfile <- args[1]

inputfile <- c("example/temp.fasta")
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
proptable <- framedata[1:512,]
UnDef <- rep(NA, 512)
proptable <- cbind(proptable, UnDef)

## Define function for handling reverse complementary pentamer
revcomp <- function(nucSeq)
  return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq)))
#rownum <- c(4, 6, 9, 14, 16, 19, 2, 3, 8, 11, 12, 13, 18, 21)
## DNAshape_retreive function2
dnashape_shapemat2 <- function(pentamer_vec) {

   rollpart2 <- c(NA)
  heltpart2 <- c(NA)
  tiltpart2 <- c(NA)
  risepart2 <- c(NA)
  shiftpart2 <- c(NA)
  slidepart2 <- c(NA)
  dim(rollpart2) <- c(1, 1)
  dim(heltpart2) <- c(1, 1)
  dim(tiltpart2) <- c(1, 1)
  dim(risepart2) <- c(1, 1)
  dim(shiftpart2) <- c(1, 1)
  dim(slidepart2) <- c(1, 1)

  ## Mapping HASh-KEY
  map <- new.env(hash = T, parent = emptyenv())
  pent5 <- rownames(framedata)
  petnew <- head(pent5,-4)
  num  <- 1:512
  for (i in seq_along(petnew))
  {
    print
    key <- petnew[i]
    map[[key]] <- num[i]
  }
  Roll_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  HelT_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  Tilt_mat <-
    array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  Rise_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  Shift_mat <- array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  Slide_mat <-
    array(numeric(), c(length(pentamer_vec) + 1, 1, 0))
  
  for (b in seq_along(pentamer_vec)) {
    checkvar <- grepl(pentamer_vec[b], rownames(proptable))
    checkvar2 <- grepl(revcomp(pentamer_vec[b]), rownames(proptable))
    if (any(checkvar)) {
      rollpart1 <- proptable[map[[pentamer_vec[b]]],][4]
      temp <- c(rollpart1[1, 1], rollpart2[1, 1])
      Roll_mat[b] <- mean(temp, na.rm = TRUE)
      rollpart2 <- proptable[map[[pentamer_vec[b]]],][5]
      
      heltpart1 <- proptable[map[[pentamer_vec[b]]],][6]
      temp <- c(heltpart1[1, 1], heltpart2[1, 1])
      HelT_mat[b] <- mean(temp, na.rm = TRUE)
      heltpart2 <- proptable[map[[pentamer_vec[b]]],][7]
      
      tiltpart1 <- proptable[map[[pentamer_vec[b]]],][9]
      temp <- c(tiltpart1[1, 1], tiltpart2[1, 1])
      Tilt_mat[b] <- mean(temp, na.rm = TRUE)
      tiltpart2 <- proptable[map[[pentamer_vec[b]]],][10]
      
      risepart1 <- proptable[map[[pentamer_vec[b]]],][14]
      temp <- c(risepart1[1, 1], risepart2[1, 1])
      Rise_mat[b] <- mean(temp, na.rm = TRUE)
      risepart2 <- proptable[map[[pentamer_vec[b]]],][15]
      
      shiftpart1 <- proptable[map[[pentamer_vec[b]]],][16]
      temp <- c(shiftpart1[1, 1], shiftpart2[1, 1])
      Shift_mat[b] <- mean(temp, na.rm = TRUE)
      shiftpart2 <- proptable[map[[pentamer_vec[b]]],][17]
      
      slidepart1 <- proptable[map[[pentamer_vec[b]]],][19]
      temp <- c(slidepart1[1, 1], slidepart2[1, 1])
      Slide_mat[b] <- mean(temp, na.rm = TRUE)
      slidepart2 <- proptable[map[[pentamer_vec[b]]],][20]
      
      
    } else if (any(checkvar2)) {
      rollpart1 <- proptable[map[[revcomp(pentamer_vec[b])]],][5]
      temp <- c(rollpart1[1, 1], rollpart2[1, 1])
      Roll_mat[b] <- mean(temp, na.rm = TRUE)
      rollpart2 <- proptable[map[[revcomp(pentamer_vec[b])]],][4]
      
      heltpart1 <- proptable[map[[revcomp(pentamer_vec[b])]],][7]
      temp <- c(heltpart1[1, 1], heltpart2[1, 1])
      HelT_mat[b] <- mean(temp, na.rm = TRUE)
      heltpart2 <- proptable[map[[revcomp(pentamer_vec[b])]],][6]
      
      tiltpart1 <- proptable[map[[revcomp(pentamer_vec[b])]],][10]
      temp <- c(tiltpart1[1, 1], tiltpart2[1, 1])
      Tilt_mat[b] <- mean(temp, na.rm = TRUE)
      tiltpart2 <- proptable[map[[revcomp(pentamer_vec[b])]],][9]
      
      risepart1 <- proptable[map[[revcomp(pentamer_vec[b])]],][15]
      temp <- c(risepart1[1, 1], risepart2[1, 1])
      Rise_mat[b] <- mean(temp, na.rm = TRUE)
      risepart2 <- proptable[map[[revcomp(pentamer_vec[b])]],][14]
      
      shiftpart1 <- proptable[map[[revcomp(pentamer_vec[b])]],][17]
      temp <- c(shiftpart1[1, 1], shiftpart2[1, 1])
      Shift_mat[b] <- mean(temp, na.rm = TRUE)
      shiftpart2 <- proptable[map[[revcomp(pentamer_vec[b])]],][16]
      
      slidepart1 <- proptable[map[[revcomp(pentamer_vec[b])]],][20]
      temp <- c(slidepart1[1, 1], slidepart2[1, 1])
      Slide_mat[b] <- mean(temp, na.rm = TRUE)
      slidepart2 <- proptable[map[[revcomp(pentamer_vec[b])]],][19]
    }else
    {
      rollpart1 <- proptable[1,][22]
      temp <- c(rollpart1[1, 1], rollpart2[1, 1])
      Roll_mat[b] <- mean(temp, na.rm = TRUE)
      
      heltpart1 <- proptable[1,][22]
      temp <- c(heltpart1[1, 1], heltpart2[1, 1])
      HelT_mat[b] <- mean(temp, na.rm = TRUE)
      
      tiltpart1 <- proptable[1,][22]
      temp <- c(tiltpart1[1, 1], tiltpart2[1, 1])
      Tilt_mat[b] <- mean(temp, na.rm = TRUE)
      
      risepart1 <- proptable[1,][22]
      temp <- c(risepart1[1, 1], risepart2[1, 1])
      Rise_mat[b] <- mean(temp, na.rm = TRUE)
      
      shiftpart1 <- proptable[1,][22]
      temp <- c(shiftpart1[1, 1], shiftpart2[1, 1])
      Shift_mat[b] <- mean(temp, na.rm = TRUE)
      
      slidepart1 <- proptable[1,][22]
      temp <- c(slidepart1[1, 1], slidepart2[1, 1])
      Slide_mat[b] <- mean(temp, na.rm = TRUE)
      
      rollpart2 <- c(NA)
      heltpart2 <- c(NA)
      tiltpart2 <- c(NA)
      risepart2 <- c(NA)
      shiftpart2 <- c(NA)
      slidepart2 <- c(NA)

    }
  }
  
  newList <-
    list(
      "Roll" = c(NA, Roll_mat, rollpart2, NA),
      "HelT" = c(NA, HelT_mat, heltpart2, NA),
      "Tilt" = c(NA, Tilt_mat, tiltpart2, NA),
      "Rise" = c(NA, Rise_mat, risepart2, NA),
      "Shift" = c(NA, Shift_mat, shiftpart2, NA),
      "Slide" = c(NA, Slide_mat, slidepart2, NA)
    )
  
  return(newList)
}

## Main script
outname <- file_path_sans_ext(basename(inputfile))

rollfilename <-
  paste(outname, ".", "Roll.mx", sep = "")
rollfilenameold <-
  paste(outname, ".", "Roll.mx.old", sep = "")
if (file.exists(rollfilename)) {
  file.rename(rollfilename, rollfilenameold)
}
heltfilename <-
  paste(outname, ".", "HelT.mx", sep = "")
heltfilenameold <-
  paste(outname, ".", "HelT.mx.old", sep = "")
if (file.exists(heltfilename)) {
  file.rename(heltfilename, heltfilenameold)
}
tiltfilename <-
  paste(outname, ".", "Tilt.mx", sep = "")
tiltfilenameold <-
  paste(outname, ".", "Tilt.mx.old", sep = "")
if (file.exists(tiltfilename)) {
  file.rename(tiltfilename, tiltfilenameold)
}
risefilename <-
  paste(outname, ".", "Rise.mx", sep = "")
risefilenameold <-
  paste(outname, ".", "Rise.mx.old", sep = "")
if (file.exists(risefilename)) {
  file.rename(risefilename, risefilenameold)
}
shiftfilename <-
  paste(outname, ".", "Shift.mx", sep = "")
shiftfilenameold <-
  paste(outname, ".", "Shift.mx.old", sep = "")
if (file.exists(shiftfilename)) {
  file.rename(shiftfilename, shiftfilenameold)
}
slidefilename <-
  paste(outname, ".", "Slide.mx", sep = "")
slidefilenameold <-
  paste(outname, ".", "Slide.mx.old", sep = "")
if (file.exists(slidefilename)) {
  file.rename(slidefilename, slidefilenameold)
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
  listshape <- dnashape_shapemat2(pentamer_vec)
  
  rollarr <-
    paste0(c(headname[1], round(as.numeric(listshape$Roll), digits = 2)), sep = "")
  write.table(
    rbind(rollarr),
    file = rollfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  heltarr <-
    paste0(c(headname[1], "NA", round(as.numeric(listshape$HelT), digits = 2), "NA"), sep = "")
  write.table(
    rbind(heltarr),
    file = heltfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  tiltarr <-
    paste0(c(headname[1], "NA", round(as.numeric(listshape$Tilt), digits = 2), "NA"), sep = "")
  write.table(
    rbind(tiltarr),
    file = tiltfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  risearr <-
    paste0(c(headname[1], "NA", round(as.numeric(listshape$Rise), digits = 2), "NA"), sep = "")
  write.table(
    rbind(risearr),
    file = risefilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  shiftarr <-
    paste0(c(headname[1], "NA", round(as.numeric(listshape$Shift), digits = 2), "NA"), sep = "")
  write.table(
    rbind(shiftarr),
    file = shiftfilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
  
  slidearr <-
    paste0(c(headname[1], "NA", round(as.numeric(listshape$Slide), digits = 2), "NA"), sep = "")
  write.table(
    rbind(slidearr),
    file = slidefilename,
    append = TRUE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    na = "NA",
    dec = ".",
    quote = FALSE
  )
}
