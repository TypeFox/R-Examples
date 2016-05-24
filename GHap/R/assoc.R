#Function: ghap.assoc
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: fit ordinary least squares for each haplotype allele

ghap.assoc<-function(
  blmm,
  haplo,
  type="HapAllele",
  gc=TRUE,
  only.active.alleles=TRUE,
  ncores=1
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
  ids <- rep(NA, times = length(blmm$residuals))
  for (i in 1:length(ids)) {
    ids[i] <- which(haplo$id == names(blmm$residuals)[i])
  }
  if (length(which(names(blmm$residuals) %in% haplo$id)) != length(blmm$residuals)) {
    stop("All ids in the GHap.blmm object must be present in the GHap.haplo object.")
  }
  
  #Check if GHap.blmm contains missing values
  if(length(na.omit(blmm$residuals)) != length(blmm$residuals)){
    stop("Missing values are not accepted.")
  }
  
  #Include weights
  blmm$w <- sqrt(1/blmm$weights)
  blmm$residuals <- blmm$w*blmm$residuals
  
  if(type == "HapAllele"){
    
    #ols iterate function
    ols.FUN <- function(j){
      x <- haplo$genotypes[j,ids]
      frq <- sum(x)/(2*length(x))
      x <- (x-mean(x))/sd(x)
      x <- blmm$w*x
      xpxi <- 1/sum(x^2)
      xpy <- sum(x*blmm$residuals)
      b <- xpxi*xpy
      se <- sqrt(var(blmm$residuals - x*b)*xpxi)
      return(c(b,se,frq))
    }
    
    #Compute haplotype regression statistics
    a <- mclapply(FUN=ols.FUN,X=which(haplo$allele.in),mc.cores = ncores)
    
    a <- data.frame(matrix(unlist(a), nrow=haplo$nalleles.in, byrow=TRUE))
    hapreg <- NULL
    hapreg$BLOCK <- haplo$block[haplo$allele.in]
    hapreg$CHR <- haplo$chr[haplo$allele.in]
    hapreg$BP1 <- haplo$bp1[haplo$allele.in]
    hapreg$BP2 <- haplo$bp2[haplo$allele.in]
    hapreg$ALLELE <- haplo$allele[haplo$allele.in]
    hapreg$BETA <- a[,1]
    hapreg$SE <- a[,2]
    hapreg$FREQ <- a[,3]
    hapreg$CHISQ.OBS <- (hapreg$BETA^2)/(hapreg$SE^2)
    hapreg$CHISQ.EXP <- qchisq(p = rank(hapreg$CHISQ.OBS)/(haplo$nalleles.in+1), df=1)
    if(gc == TRUE){
      dev <- sd(hapreg$CHISQ.OBS)
      samp <- which(hapreg$CHISQ.OBS < 3*dev)
      lambda <- lm(hapreg$CHISQ.OBS[samp] ~ hapreg$CHISQ.EXP[samp])
      lambda <- lambda$coefficients[2]
      hapreg$CHISQ.OBS <- hapreg$CHISQ.OBS/lambda
    }
    hapreg$logP <- -1*pchisq(q = hapreg$CHISQ.OBS, df = 1, lower.tail=FALSE, log.p = TRUE)/log(10)
    hapreg <- data.frame(hapreg,stringsAsFactors = FALSE)
  } else if(type == "HapBlock") {
    
    #ols iterate function for HapBlock
    ols.FUN <- function(j){
      hapalleles <- which(haplo$block == j)
      nalleles <- length(hapalleles)
      X <- as.matrix(haplo$genotypes[hapalleles,ids])
      if(nalleles > 1){
        X <- t(X)
        X <- scale(X)
        X <- blmm$w*X
        Xp <- t(X)
      }else{
        X <- scale(X)
        X <- blmm$w*X
        Xp <- t(X)
      }
      XpXi <- solve(Xp%*%X)
      Xpy <- Xp%*%blmm$residuals
      b <- XpXi%*%Xpy
      Xb <- X%*%b
      RSS <- sum((blmm$residuals - Xb)^2)
      ESS <- sum((Xb - mean(blmm$residuals))^2)
      nu1 <- nalleles
      nu2 <- length(blmm$residuals)-nalleles
      ftest <- (ESS/nu1)/(RSS/nu2)
      return(c(nalleles,nu1,nu2,ftest))
    }
    
    #Get unique blocks
    blocks <- unique(data.frame(haplo$block,haplo$chr,haplo$bp1,haplo$bp2,stringsAsFactors = FALSE))
    colnames(blocks) <- c("BLOCK","CHR","BP1","BP2")
    
    #Compute haplotype regression statistics
    a <- mclapply(FUN=ols.FUN,X=blocks$BLOCK,mc.cores = ncores)
    a <- data.frame(matrix(unlist(a), nrow=nrow(blocks), byrow=TRUE))
    
    hapreg <- NULL
    hapreg$BLOCK <- blocks$BLOCK
    hapreg$CHR <- blocks$CHR
    hapreg$BP1 <- blocks$BP1
    hapreg$BP2 <- blocks$BP1
    hapreg$N.ALLELES <- as.integer(a[,1])
    hapreg$F.TEST <- a[,4]
    hapreg$logP <- -1*pf(q = hapreg$F.TEST, df1 = a[,2], df2 = a[,3], lower.tail=FALSE, log.p = TRUE)/log(10)
    hapreg <- data.frame(hapreg,stringsAsFactors = FALSE)
  } else {
    stop("Argument type must be 'HapAllele' or 'HapBlock'")
  }
  
  
  #Return object
  return(hapreg)
  
}