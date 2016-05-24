#Function: ghap.kinship
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Compute haplotype covariance matrix

ghap.kinship<-function(
  haplo,
  weights=NULL,
  batchsize=500,
  only.active.samples=TRUE,
  only.active.alleles=TRUE,
  verbose=TRUE
){
  
  #Check if haplo is a GHap.haplo object
  if(class(haplo) != "GHap.haplo"){
    stop("Argument haplo must be a GHap.haplo object.")
  }
  
  #Check if inactive alleles and samples should be reactived
  if(only.active.alleles == FALSE){
    haplo$allele.in <- rep(TRUE,times=haplo$nalleles)
    haplo$nalleles.in<-length(which(haplo$allele.in))
  }
  if(only.active.samples == FALSE){
    haplo$id.in <- rep(TRUE,times=haplo$nsamples)
    haplo$nsamples.in<-length(which(haplo$id.in))
  }
  
  #Check weights
  if (is.null(weights) == TRUE) {
    weights <- rep(1,times=haplo$nalleles.in)
  }
  if (length(weights) != haplo$nalleles.in) {
    stop("Vector of weights must have the same length as the number of haplotype alleles.")
  }
  
  #Generate batch index
  id1<-seq(1,haplo$nalleles,by=batchsize)
  id2<-(id1+batchsize)-1
  id1<-id1[id2<=haplo$nalleles]
  id2<-id2[id2<=haplo$nalleles]
  id1 <- c(id1,id2[length(id2)]+1)
  id2 <- c(id2,haplo$nalleles)
  if(id1[length(id1)] > haplo$nalleles){
    id1 <- id1[-length(id1)]; id2 <- id2[-length(id2)]
  }
  
  #Log message
  if(verbose == TRUE){
    cat("Processing ", haplo$nalleles, " haplotype alleles in:\n", sep="")
    batch <- table((id2-id1)+1)
    for(i in 1:length(batch)){
      cat(batch[i]," batches of ",names(batch[i]),"\n",sep="")
    }
    cat("Inactive alleles will be ignored.\n")
  }
  
  #Initialize kinship matrix
  cat("Preparing",haplo$nsamples.in,"x",haplo$nsamples.in,"kinship matrix.\n")
  K <- matrix(data = 0, nrow = haplo$nsamples.in, ncol = haplo$nsamples.in)
  
  #Kinship iterate function
  activealleles <- which(haplo$allele.in)
  kinship.FUN<-function(i){
    slice <- id1[i]:id2[i]
    slice <- slice[slice %in% activealleles]
    Ztmp <- haplo$genotypes[slice,haplo$id.in]
    Ztmp.mean <- apply(X = Ztmp,MARGIN = 1,FUN = mean)
    #Ztmp.sd <- apply(X = Ztmp,MARGIN = 1,FUN = sd)
    #Ztmp <- (Ztmp - Ztmp.mean)/Ztmp.sd
    Ztmp <- (Ztmp - Ztmp.mean)
    K <- K + t(Ztmp)%*%(weights[slice]*Ztmp)
    return(K)
  }
  
  #Loop by batch
  for(i in 1:length(id1)){
    K <- kinship.FUN(i)
    if(verbose == TRUE){
      cat(id2[i], "haplotype alleles processed.\r")
    }
  }
  
  #Scale kinship matrix
  #K <- K/haplo$nalleles.in
  varfun <- function(j) return(var(haplo$genotypes[j,]))
  q <- sum(unlist(mclapply(X=which(haplo$allele.in),FUN = varfun)))
  K <- K/q
  colnames(K) <- haplo$id[haplo$id.in]
  rownames(K) <- colnames(K)
  
  #Return output
  return(K)
  
}