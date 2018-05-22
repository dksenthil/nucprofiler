## Load Librarys
library(data.table)
library(seqinr)
library(stringi)
library("Biostrings")
library(tools)
library(MASS)
library(rtracklayer)
library(GenomeInfoDb)



## get commandline options as input
#args <- c("SingleSeqsample.fa", "ProT,MGW,HelT1") # comment thi line if running in commandline
args <-
  c("GCGTA.fasta","Roll,MGW,EP,ProT,Tilt,Roll,HelT,Shift,Slide,Rise,Shear,Stretch,Stagger,Buckle,PropT,Opening") # comment thi line if running in commandline
    #args <- c("GCGTA.fasta","MGW,EP,Tilt,Roll,HelT,Shift,Slide,Rise,Shear,Stretch,Stagger,Buckle,PropT,Opening"
    #args <- c("GCGTA.fasta", "MGW") # comment thi line if running in commandline
    #args <- c("outchr1.fa", "Roll") # comment thi line if running in commandline
    
    #args <- commandArgs(trailingOnly = TRUE)
    # test if there is at least one argument: if not, return an error
    if (length(args) < 2)
      stop ("Usage: Rscipt  run_DNAshapeProfiler.R <Fastafile> <CommaSeparated Property list>")
    
    ## Read property table
    mydata <- fread("../data/Table_property_from_DiProDB.csv")
    mytable <-
      mydata[, list(PropertyName,
                    AA,
                    AC,
                    AG,
                    AT,
                    CA,
                    CC,
                    CG,
                    CT,
                    GA,
                    GC,
                    GG,
                    GT,
                    TA,
                    TC,
                    TG,
                    TT)]
    rohdata <- fread("../data/DNA_shape_query_table.csv")
    proptable <-
      rohdata[, list(
        Pentamer,
        MGW,
        ProT,
        Roll1,
        Roll2,
        HelT1,
        HelT2,
        Stretch,
        Tilt1,
        Tilt2,
        Buckle,
        Shear,
        Opening,
        Rise1,
        Rise2,
        Shift1,
        Shift2,
        Stagger,
        Slide1,
        Slide2,
        EP
      )]
    
    ##  Load the DNAshape short sequence
    fastafilename <- readDNAStringSet(args[1], format = "fasta")
 #   outname <- basename(args[1])
  #  print(outname)
    ## DNA table pentamer calculation
    
    pentamer <- vector()
    sequence = paste(fastafilename)
    seqvec <- strsplit(sequence, "")
    for (i in 1:(nchar(sequence) - 4)) {
      pentamer[i] <- paste0(seqvec[[1]][i:(i + 4)], collapse = "")
    }
    
    ## Get the table property in the loop:
    ## Assumes DNAshape property name strict and only 20 property names with suffix (1 or 2)
    
    proplist <-
      c(
        "Name",
        "MGW",
        "ProT",
        "Stretch",
        "Buckle",
        "Shear",
        "Opening",
        "Stagger",
        "Shift",
        "Rise",
        "Tilt",
        "Slide",
        "Roll",
        "HelT",
        "EP"
      )
    
    getprop <- unlist(strsplit(args[2], ",", fixed = TRUE))
    revcomp <- function(nucSeq)
      # uses library stringi
      return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq)))
    
    outmat<-vector()  
    neighborval <- c(NA)
    for (q in 1:length(getprop)) {
      if (getprop[q] == "Roll" ||
          getprop[q] ==  "HelT" ||
          getprop[q] == "Tilt" ||
          getprop[q] == "Rise" ||
          getprop[q] == "Shift" || getprop[q] == "Slide") {
        # assign the column to pick
        if (getprop[q] == "Roll") {
          rownum <- 4
          n <- 1
        }
        if (getprop[q] == "HelT") {
          rownum <- 6
          n <- 1
        }
        if (getprop[q] == "Tilt") {
          rownum <- 9
          n <- 1
        }
        if (getprop[q] == "Rise") {
          rownum <- 14
          n <- 1
        }
        if (getprop[q] == "Shift") {
          rownum <- 16
          n <- 1
        }
        if (getprop[q] == "Slide") {
          rownum <- 19
          n <- 1
        }
        ##
        
        for (i in 1:length(pentamer)) {
          #    for(i in 1:3){
          for (j in 1:512) {
            if (pentamer[i] == proptable[[1]][j]) {
              if (i == 1) {
                outmat[i] <- c(proptable[[rownum]][j])
                neighborval <- c(proptable[[rownum + 1]][j])
                print("first", neighborval)
                
              } else if (i > 1 && i < length(pentamer)) {
                if(is.na(neighborval)){
                  outmat[i] <- c(proptable[[rownum]][j])
                  neighborval <- proptable[[rownum + 1]][j]
                  print ("ok1")
                }else {
                avgval <- mean(c(neighborval, proptable[[rownum]][j]))
                outmat[i] <- c(avgval)
                neighborval <- proptable[[rownum + 1]][j]
                print ("ok2")
                }
              } else if (i == length(pentamer))
              {
                if(is.na(neighborval)){
                  outmat[i] <- c(proptable[[rownum]][j])
                  outmat[i + 1]  <- c(proptable[[rownum + 1]][j])
                }else {
                avgval <- mean(c(neighborval, proptable[[rownum]][j]))
                outmat[i] <- c(avgval)
                outmat[i + 1]  <- c(proptable[[rownum + 1]][j])
                }
              }
            } else  if (revcomp(pentamer[i]) == proptable[[1]][j]) {
              if (i == 1) {
                outmat[i] <- c(proptable[[rownum + 1]][j])
                neighborval <- proptable[[rownum]][j]
                print(neighborval)
                print ("ok")
              } else if (i > 1 && i < length(pentamer)) {
                if(is.na(neighborval)){
                  outmat[i] <- c(proptable[[rownum+1]][j])
                }else {              
                avgval <- mean(c(neighborval, proptable[[rownum + 1]][j]))
                outmat[i] <- c(avgval)
                neighborval <- proptable[[rownum]][j]
                print ("ok")
                }
              } else if (i == length(pentamer))
              {
                if(is.na(neighborval)){
                  outmat[i] <- c(proptable[[rownum+1]][j])
                  outmat[i + 1]  <- c(proptable[[rownum]][j])
                  print("here1")
                }else { 
                avgval <- mean(c(neighborval, proptable[[rownum + 1]][j]))
                outmat[i] <- c(avgval)
                outmat[i + 1]  <- c(proptable[[rownum]][j])
                print("here2")
                }
              }
              
            }else if (grepl("N", pentamer[i])){
   #           if(i < length(pentamer)){
                outmat[i] <- c(neighborval)
                neighborval <- c(NA)
                print("Yep")
                
 #             }if (i == length(pentamer)){
  #            neighborval <- c(NA)   
 #             avgval <- mean(c(neighborval, proptable[[rownum + 1]][j]))
              #outmat[i]  <- c(neighborval)
            }
            }

          }
          
          # Shift, Rise, Tilt, Slide, Roll and HelT
          # Write output matrix in  a WIG file format to a file.
          newvec <- vector()
          outname <- file_path_sans_ext(basename(args[1]))
          for (i in 1:length(outmat))
          {
            newvec[i] <- toString(round((i + n), digits = 0))
          }
          
          outmat2 <-
            paste(newvec, round(outmat, digits = 2), sep = " ") # space separted columns in matrix
          print(outmat2)
          chromo <- strsplit(names(fastafilename), " ")
          chromosomename <- chromo[[1]][1]
          chrom <- paste0("variableStep chrom=",chromosomename)
          write(chrom, file = paste(outname, ".", getprop[q],".wig", sep = ""), sep = "")         
          write(outmat2,
                       file = paste(outname, ".", getprop[q],".wig", sep = ""),
                       sep = "", append = TRUE)
          
        }
    
   ##----------------use subroutine-------
   else if (getprop[q] == "MGW" ||
            getprop[q] == "ProT" ||
            getprop[q] == "Stretch" ||
            getprop[q] == "Shear" ||
            getprop[q] == "Buckle" ||
            getprop[q] == "Opening" ||
            getprop[q] == "Stagger" || getprop[q] == "EP") {
     ##
     if (getprop[q] == "MGW") {
       rownum <- 2
       n <- 2
     }
     if (getprop[q] == "ProT") {
       rownum <- 3
       n <- 2
     }
     if (getprop[q] == "Stretch") {
       rownum <- 8
       n <- 2
     }
     if (getprop[q] == "Shear") {
       rownum <- 12
       n <- 2
     }
     if (getprop[q] == "Buckle") {
       rownum <- 11
       n <- 2
     }
     if (getprop[q] == "Opening") {
       rownum <- 13
       n <- 2
     }
     if (getprop[q] == "Stagger") {
       rownum <- 18
       n <- 2
     }
     if (getprop[q] == "EP") {
       rownum <- 21
       n <- 2
     }
     
     # run the property value pick loop
     outmat <- vector()
     neighborval <- vector()
     length(outmat) <- length(pentamer)
     length.outmat = length(pentamer) + 1
     # neighborval <- c(0)
     
     ## String reverse and check missing values
     revcomp <- function(nucSeq)
    uses library stringi
       return(stri_reverse(chartr("acgtACGT", "tgcaTGCA", nucSeq)))
     
     
     for (i in 1:length(pentamer)) {
       for(i in 1:3){
       for (j in 1:512) {
         if (pentamer[i] == proptable[[1]][j])
         {
           outmat[i]  <- c(proptable[[rownum]][j])
         } else  if (revcomp(pentamer[i]) == proptable[[1]][j]) {
           outmat[i]  <- c(proptable[[rownum]][j])
         }
         
       }
     }
     
     # Shift, Rise, Tilt, Slide, Roll and HelT
     # Write output matrix in  a WIG file format to a file.
     newvec <- vector()
     outname <- file_path_sans_ext(basename(args[1]))
     for (i in 1:length(outmat))
     {
       newvec[i] <- toString(round((i + n), digits = 0))
     }
     
     outmat2 <-
       paste(newvec, round(outmat, digits = 2), sep = " ") # space separted columns in matrix
     print(outmat2)
     chromo <- strsplit(names(fastafilename), " ")
     chromosomename <- chromo[[1]][1]
     chrom <- paste0("variableStep chrom=",chromosomename)
     write(chrom, file = paste(outname, ".", getprop[q],".wig", sep = ""), sep = "")
     write(outmat2, file = paste(outname, ".", getprop[q],".wig", sep = ""), sep = "", append = TRUE)
     
   }
    }

# convert WIG to BIGWIG
    
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
      
    
    
    
