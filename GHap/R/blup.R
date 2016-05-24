#Function: ghap.gblup
#License: GPLv3 or later
#Modification date: 6 Apr 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: calculate GBLUP solution for each haplotype allele

ghap.blup<-function(
  blmm, 
  haplo,
  weights = NULL,
  nperm = 1,
  only.active.alleles = TRUE, 
  ncores = 1
){
  
  #General data check
  if (class(blmm) != "GHap.blmm") {
    stop("Argument blmm must be a GHap.blmm object.")
  }
  if (class(haplo) != "GHap.haplo") {
    stop("Argument haplo must be a GHap.haplo object.")
  }
  if (only.active.alleles == FALSE) {
    haplo$allele.in <- rep(TRUE, times = haplo$nalleles)
    haplo$nalleles.in <- length(which(haplo$allele.in))
  }
  if (is.null(weights) == TRUE) {
    weights <- rep(1,times=haplo$nalleles.in)
  }
  if (length(weights) != haplo$nalleles.in) {
    stop("Vector of weights must have the same length as the number of haplotype alleles.")
  }
  ids <- rep(NA, times = length(blmm$k))
  for (i in 1:length(ids)) {
    ids[i] <- which(haplo$id == names(blmm$k)[i])
  }
  if (length(which(names(blmm$k) %in% haplo$id)) != length(blmm$k)) {
    stop("All ids in the GHap.blmm object must be present in the GHap.haplo object.")
  }
  
  #Compute variance
  varfun <- function(j) return(var(haplo$genotypes[j,]))
  sumvar <- sum(unlist(mclapply(X=which(haplo$allele.in),FUN = varfun,  mc.cores = ncores)))
  
  #Compute frequencies
  freqfun <- function(j) return(sum(haplo$genotypes[j,])/(2*haplo$nsamples.in))
  freq <- unlist(mclapply(X=which(haplo$allele.in),FUN = freqfun,  mc.cores = ncores))
  
  #Main BLUP function
  gblup.FUN <- function(j) {
    x <- haplo$genotypes[j, ids]
    x <- x - mean(x)
    b <- sum(weights[j]*x*blmm$k)
    b <- b/sumvar
    varxb <- var(x*b)
    return(c(b,varxb))
  }
  
  #Compute effects
  a <- mclapply(FUN = gblup.FUN, X = which(haplo$allele.in), mc.cores = ncores)
  a <- data.frame(matrix(unlist(a), nrow=haplo$nalleles.in, byrow=TRUE))
  
  #Output data
  hapreg <- NULL
  hapreg$BLOCK <- haplo$block[haplo$allele.in]
  hapreg$CHR <- haplo$chr[haplo$allele.in]
  hapreg$BP1 <- haplo$bp1[haplo$allele.in]
  hapreg$BP2 <- haplo$bp2[haplo$allele.in]
  hapreg$ALLELE <- haplo$allele[haplo$allele.in]
  hapreg$SCORE <- a[,1]
  hapreg$FREQ <- freq
  hapreg$VAR <- a[,2]
  hapreg$pVAR <- hapreg$VAR/sum(hapreg$VAR)
  
  #Permutation test (optional)
  if(nperm > 1){
    cat("A permutation procedure with",nperm,"randomizations will be performed.\n")
    #Re-define BLUP function to randomize effects
    #Main BLUP function
    gblup.FUN <- function(j) {
      x <- haplo$genotypes[j, ids]
      x <- x - mean(x)
      b <- sum(weights[j]*x*blmm$k)
      b <- b/sumvar
      return(b)
    }
    hapreg$P <- 0
    #Permutation iteration
    for(i in 1:nperm){
      cat("Permutation number:",i,"\r")
      blmm$k <- sample(blmm$k,size=length(blmm$k),replace=FALSE)
      a <- mclapply(FUN = gblup.FUN, X = which(haplo$allele.in), mc.cores = ncores)
      a <- unlist(a)
      a <- as.numeric(abs(max(a)) > abs(hapreg$SCORE))
      hapreg$P <- hapreg$P + a
    }
    hapreg$P <- hapreg$P/nperm
    hapreg$P[hapreg$P == 0] <- 1/nperm
  }
  
  #Return results
  hapreg <- data.frame(hapreg, stringsAsFactors = FALSE)
  return(hapreg)
}

