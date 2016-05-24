# this script takes the Beagle output with the 1 2 Allele coding, imports the Variant Master Key and substitutes then the 
# Alles with the corrsepsonding nucleotide.

subPhased <- function(file=NULL, vmmk = NULL, out=NULL, chunkSize=100000, verbose=TRUE, removeInsertions=TRUE){

  # Basic input checks
    if(is.null(file)) stop("No *.phased file given.")
    if(is.null(vmmk)) stop("No Variant-Map-Master-Key file given.")
    if(is.null(out)) out <- paste("Allele-",file,sep="")
  # Avoid E-Notation 
    options("scipen"=100, "digits"=4)
  # Read in the variant master key:
    master <- read.table(vmmk, stringsAsFactors=FALSE)
    if(verbose) cat("I read the master key file",date(),"\n")
  # Determine the number of repetitions
  # nrows <- countLines(file) - 1
  nrows <- system(paste("wc -l",file), intern=TRUE) 
  nrows <- as.numeric(strsplit(nrows, " ")[[1]][1]) - 1
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
      if(verbose) cat("I read the phased chunk",chunkRun,"/",nchunks,":",date(),"\n")
    # Check if the order is correct, therwise stop the script!
      if(sum(phased$id==master$V1[((chunkRun-1)*chunkSize+1):((chunkRun-1)*chunkSize+readRows)])!=nrow(phased)) stop("The order doesn't match in the files.")
      phased12 <- phased[,1:2]
      phased <- phased[,-c(1,2)]
      phased <- cbind(master[((chunkRun-1)*chunkSize+1):((chunkRun-1)*chunkSize+readRows),3:4],phased)
    # Now change the dose table
  #    cat("I prepared everything",date(),"\n")
      phased <- t(apply(phased,1,function(x){c <- x[1:2]; y <- x[-(1:2)]; ifelse(y==1,c[1],c[2])}))
  #    cat("I run the apply",date(),"\n")
      phased <- cbind(phased12,phased)
  #    cat("Right format",date(),"\n")
    # Remove the Insertions
      if(removeInsertions) phased <- phased[which(master$V7[((chunkRun-1)*chunkSize+1):((chunkRun-1)*chunkSize+readRows)]=="S"),]
  #    cat("Removed I entries",date(),"\n")
    # write out the new Dose file
      if(chunkRun==1){
        phasedCN[1:2] <- colnames(phased)[1:2]
        colnames(phased) <- phasedCN
        write.table(phased,out, col.names=TRUE, row.names=FALSE, quote=FALSE)
      } else {
        write.table(phased,out, col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
      }   
    }
}