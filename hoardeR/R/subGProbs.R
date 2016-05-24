# this script takes the Beagle output with the 1 2 Allele coding, imports the Variant Master Key and substitutes then the 
# Alles with the corrsepsonding nucleotide.

subGprobs <- function(file=NULL, vmmk=NULL, out=NULL, chunkSize=100000, removeInsertions=TRUE, verbose = TRUE, writeOut=TRUE){

  # Basic input checks
    if(is.null(file)) stop("No *.gprobs file given.")
    if(is.null(vmmk)) stop("No Variant-Map-Master-Key file given.")
    if(is.null(out)) out <- paste("Allele-",file,sep="")

  # Avoid E-Notation 
    options("scipen"=100, "digits"=4)
  # Read in the variant master key:
    master <- read.table(vmmk, stringsAsFactors=FALSE)
    cat("I read the master key file",date(),"\n")
  # Determine the number of repetitions
    # nrows <- countLines(file) - 1
    nrows <- system(paste("wc -l",file), intern=TRUE) 
    nrows <- as.numeric(strsplit(nrows, " ")[[1]][1]) - 1
    if(verbose) cat("Number of rows:", nrows,"\n")
    nchunks <- ceiling(nrows / chunkSize)

  # Now process chunk by chunk
    for(chunkRun in 1:nchunks){
      gc()
      skipThis <- chunkSize*(chunkRun-1)
      ifelse((nrows - skipThis) > chunkSize, readRows <- chunkSize, readRows <- nrows-skipThis)
      if(chunkRun==1){
        phased <- read.table(file, check.names=FALSE, header=TRUE, stringsAsFactors=FALSE, skip=skipThis, nrows=readRows)
        phasedCN <- colnames(phased)
      } else {
        phased <- read.table(file, header=FALSE, stringsAsFactors=FALSE, skip=skipThis+1, nrows=readRows)
        colnames(phased) <- phasedCN
      }
      cat("I read the GProbs chunk",chunkRun,"/",nchunks,":",date(),"\n")
    # Check if the order is correct, therwise stop the script!
      if(sum(phased$marker==master$V1[((chunkRun-1)*chunkSize+1):((chunkRun-1)*chunkSize+readRows)])!=nrow(phased)) stop("The order doesn't match in run",chunkRun)
    # substitute the columns
      phased$alleleA <- master$V3[((chunkRun-1)*chunkSize+1):((chunkRun-1)*chunkSize+readRows)]
      phased$alleleB <- master$V4[((chunkRun-1)*chunkSize+1):((chunkRun-1)*chunkSize+readRows)]
      
    # Remove the Insertions
      if(removeInsertions) phased <- phased[which(master$V7[((chunkRun-1)*chunkSize+1):((chunkRun-1)*chunkSize+readRows)]=="S"),]

    # write out the new Dose file
    if(writeOut){
      if(chunkRun==1){
        write.table(phased,out, col.names=TRUE, row.names=FALSE, quote=FALSE)
      } else {
        write.table(phased,out, col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
      } 
    }
    
    }  
}

