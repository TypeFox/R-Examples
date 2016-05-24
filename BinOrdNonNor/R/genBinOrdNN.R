genBinOrdNN <-
function(n, plist, mean.vec, var.vec, skew.vec, kurto.vec, no.bin, no.ord, no.NN, cmat.star) 
{ 
  
  if (missing(n)){
    stop("n was not specified! \n")
  }
  
  if (missing(cmat.star)) {
    stop("The intermediate correlation matrix was not specified! \n")    
  }
  
  if (no.bin > 0){
    for (i in 1:no.bin){
      if (length(plist[[i]])>1) {
        warning("The probability vector for the ", i, "-th binary variable contains more than one value!\n") 
      }
    }    
  }  
  
  if (no.ord > 0){
    for (i in 1:no.ord){
      if (length(plist[[no.bin+i]])<2) {
        warning("The probability vector for the ", eval(no.bin+i), "-th ordinal variable contains only one value!\n")
      }
    }    
  }
  
  no.binord <- no.bin + no.ord
    
  if (ncol(cmat.star) != (no.binord + no.NN)){
    stop("Dimension of intermediate correlation matrix cmat.star des not match the number of variables!\n")
  }  
  
  if (no.binord > 0) {
    if (length(plist) != no.binord) {
      stop("Dimension of the probability vector does not match the number of binary and ordinal variables!\n")
    }  
  }
  
  if (no.NN > 0){
    if (length(skew.vec) != no.NN) {
      stop("Length of the skewness vector does not match the number of non-normal variables!\n")
    }
    
    if (length(kurto.vec) != no.NN) {
      stop("Length of the kurtosis vector does not match the number of non-normal variables!\n")
    }
  
    
    if (length(mean.vec) != no.NN) {
      stop("Length of the mean vector does not match the number of continuous variables!\n")
    }
    
    if (length(var.vec) != no.NN) {
      stop("Length of the variance vector does not match the number of continuous variables!\n")
    }
  }
    
  if (no.NN == 0) {
    YY <- ordsample(n, plist, cmat.star, cormat = "continuous")
  }
  
  if (no.binord == 0) {
    coef <- Fleishman.coef.NN(skew.vec,kurto.vec)
    XX <- rmvnorm(n, rep(0, ncol(cmat.star)), cmat.star)
    YY <- NULL
    for (i in 1:no.NN){
      X <- cbind(1,XX[,i], XX[,i]^2, XX[,i]^3)
      Y <- X%*%coef[i,]*sqrt(var.vec[i-length(plist)])+mean.vec[i-length(plist)]
      YY <- cbind(YY,Y)          
    }
  }
  
  if (no.NN > 0 & no.binord > 0) {
    XX <- rmvnorm(n, rep(0, ncol(cmat.star)), cmat.star)
    YY <- NULL
    
    for (i in 1:no.binord) {
      OO <- ordinalize(plist[[i]], XX[, i])
      YY <- cbind(YY, OO)
      rm(OO)
    }
    
    coef <- Fleishman.coef.NN(skew.vec,kurto.vec)
    
    for (i in (no.binord + 1):(no.binord + no.NN)) {
      X <- cbind(1,XX[,i], XX[,i]^2, XX[,i]^3)
      Y <- X%*%coef[i-length(plist),]*sqrt(var.vec[i-length(plist)])+mean.vec[i-length(plist)]
      YY <- cbind(YY,Y)          
    }
    
    rm(XX)
  }
  
  colnames(YY) <- NULL
  return(YY)
}
