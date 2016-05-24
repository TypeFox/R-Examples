#Function: ghap.mergehaplo
#License: GPLv3 or later
#Modification date: 15 Apr 2016
#Written by: Marco Milanesi & Yuri Tani Utsunomiya
#Contact: marco.milanesi.mm@gmail.com, ytutsunomiya@gmail.com
#Description: Merge GHap.haplo objects 

ghap.mergehaplo<-function(
  haplo.1,  haplo.2,  type, only.active.markers=TRUE, only.active.samples=TRUE, verbose=TRUE
){
  
  #Check if haplo.1 is a GHap.haplo object
  if(class(haplo.1) != "GHap.haplo"){
    stop("Argument haplo.1 must be a GHap.haplo object.")
  }
  #Check if haplo.2 is a GHap.haplo object
  if(class(haplo.2) != "GHap.haplo"){
    stop("Argument haplo.2 must be a GHap.haplo object.")
  }
  
  #Check if inactive markers and samples should be reactived
  if(only.active.markers == FALSE){
    haplo.1$allele.in <- rep(TRUE,times=haplo.1$nalleles)
    haplo.2$allele.in <- rep(TRUE,times=haplo.2$nalleles)
    haplo.1$nalleles.in <- haplo.1$nalleles
    haplo.2$nalleles.in <- haplo.2$nalleles
  }
  if(only.active.samples == FALSE){
    haplo.1$id.in <- rep(TRUE,times=haplo.1$nsamples)
    haplo.2$id.in <- rep(TRUE,times=haplo.2$nsamples)
    haplo.1$nsamples.in <- haplo.1$nsamples
    haplo.2$nsamples.in <- haplo.2$nsamples
  }
  
  # create temporary block_allele variable
  tmpHB1 <- paste(haplo.1$block[haplo.1$allele.in], haplo.1$allele[haplo.1$allele.in], sep="")
  tmpHB2 <- paste(haplo.2$block[haplo.2$allele.in], haplo.2$allele[haplo.2$allele.in], sep="")
  
  # Start merging
  if (type == "HapAllele") {
    
    #Check if there is duplicate HapAllele
    if(length(which(tmpHB1 %in% tmpHB2)) != 0 | length(which(tmpHB2 %in% tmpHB1)) != 0){
      stop("Duplicated HapAlleles are not allowed.")
    }
    
    #Check if IDs are the same and in the same order
    if(identical(haplo.1$id[haplo.1$id.in],haplo.2$id[haplo.2$id.in]) == FALSE){
      stop("IDs in the two GHap.phase objects are not the same or are not in the same order.")
    }
    
    #Create the merged GHap.haplo object
    if(verbose == TRUE){
      cat("\nCreating the new GHap.haplo object... ")
    }
    haplo<-NULL
    haplo$nsamples <- haplo.1$nsamples.in
    haplo$nalleles <- (haplo.1$nalleles.in+haplo.2$nalleles.in)
    haplo$nsamples.in <- haplo.1$nsamples.in
    haplo$nalleles.in <- (haplo.1$nalleles.in+haplo.2$nalleles.in)
    haplo$pop <- haplo.1$pop[haplo.1$id.in]
    haplo$id <- haplo.1$id[haplo.1$id.in]
    haplo$id.in<-rep(TRUE,haplo$nsamples.in)
    haplo$chr<-c(haplo.1$chr[haplo.1$allele.in],haplo.2$chr[haplo.2$allele.in])
    haplo$block<-c(haplo.1$block[haplo.1$allele.in],haplo.2$block[haplo.2$allele.in])
    haplo$bp1<-c(haplo.1$bp1[haplo.1$allele.in],haplo.2$bp1[haplo.2$allele.in])
    haplo$bp2<-c(haplo.1$bp2[haplo.1$allele.in],haplo.2$bp2[haplo.2$allele.in])
    haplo$allele<-c(haplo.1$allele[haplo.1$allele.in],haplo.2$allele[haplo.2$allele.in])
    haplo$allele.in<-rep(TRUE,haplo$nalleles.in)
    tmp1 <- as.matrix(haplo.1$genotypes[haplo.1$allele.in,haplo.1$id.in]) # export to a matrix first haplo data
    tmp2 <- as.matrix(haplo.2$genotypes[haplo.2$allele.in,haplo.2$id.in]) # export to a matrix second haplo data
    options(bigmemory.typecast.warning=FALSE)
    haplo$genotypes <- as.big.matrix(rbind(tmp1,tmp2), type="char")
    rm(tmp1,tmp2)
    if(verbose == TRUE){
      cat("Done.\n")
    }
    
  } else if (type == "individual") {
    
    #Check if there is duplicate IDs
    if(length(which(haplo.1$id[haplo.1$id.in] %in% haplo.2$id[haplo.2$id.in])) != 0){
      stop("Duplicated IDs are not allowed.")
    }
    
    #Check if HapAlelles are the same and in the same order
    if(identical(tmpHB1,tmpHB2) == FALSE){
      stop("HapAlleles in the two GHap.phase objects are not the same or are not in the same order.")
    }
    
    #Create the merged GHap.haplo object
    if(verbose == TRUE){
      cat("\nCreating the new GHap.haplo object... ")
    }
    haplo<-NULL
    haplo$nsamples <- (haplo.1$nsamples.in+haplo.2$nsamples.in)
    haplo$nalleles <- haplo.1$nalleles.in
    haplo$nsamples.in <- (haplo.1$nsamples.in+haplo.2$nsamples.in)
    haplo$nalleles.in <- haplo.1$nalleles.in
    haplo$pop <- c(haplo.1$pop[haplo.1$id.in],haplo.2$pop[haplo.2$id.in])
    haplo$id <- c(haplo.1$id[haplo.1$id.in],haplo.2$id[haplo.2$id.in])
    haplo$id.in<-rep(TRUE,haplo$nsamples.in)
    haplo$chr<-haplo.1$chr[haplo.1$allele.in]
    haplo$block<-haplo.1$block[haplo.1$allele.in]
    haplo$bp1<-haplo.1$bp1[haplo.1$allele.in]
    haplo$bp2<-haplo.1$bp2[haplo.1$allele.in]
    haplo$allele<-haplo.1$allele[haplo.1$allele.in]
    haplo$allele.in<-rep(TRUE,haplo$nalleles.in)
    tmp1 <- as.matrix(haplo.1$genotypes[haplo.1$allele.in,haplo.1$id.in]) # export to a matrix first haplo data
    tmp2 <- as.matrix(haplo.2$genotypes[haplo.2$allele.in,haplo.2$id.in]) # export to a matrix second haplo data
    options(bigmemory.typecast.warning=FALSE)
    haplo$genotypes <- as.big.matrix(cbind(tmp1,tmp2), type="char")
    rm(tmp1,tmp2)
    if(verbose == TRUE){
      cat("Done.\n")
    }
    
  }else{
    stop("Argument type must be 'individual' or 'HapAllele'")
  }
  
  
  #Return ghap object
  class(haplo) <- "GHap.haplo"
  if(verbose == TRUE){
    cat("Your GHap.haplo object was successfully merged without apparent errors.\n\n")
  }
  return(haplo)
}
