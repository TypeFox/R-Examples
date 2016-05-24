ExpandCoefBasis<-function(coef, ncol, splinebasis, expand=rep(TRUE, ncol), value=0){
# splinebasis : object of class "BSplineBasis",, "MSplineBasis",, "TPSplineBasis",
# transform a vector of coef into a matrix of coef for spline evluation
# if expand[i] ==FALSE, coef contains the nTbasis spline coef of the var i
# if expand[i] ==True, coef contains the nTbasis-1 spline coef of the var i, thus add "value"  in the first col
  stopifnot(is.integer(ncol))
  if(length(value)>1){
    stopifnot(length(value)==ncol)
  }
  else {
    value <- rep(value, ncol)
  }
  mcoef <- matrix(value, ncol=ncol, nrow=splinebasis@nbases + splinebasis@log, byrow=TRUE)
  if(length(expand)>1){
    stopifnot(length(expand)==ncol)
    first <- 1 + expand
    nb <- splinebasis@nbases + splinebasis@log - expand 
  }
  else {
    first <- 1 + rep(expand, ncol)
    nb <- splinebasis@nbases + splinebasis@log - rep(expand, ncol)
  }  
  beg <- 0
  for(i in 1:ncol){
    end <- beg + nb[i]
    beg <- beg+1
    mcoef[first[i]:(splinebasis@nbases + splinebasis@log),i] <- coef[beg:end]
    beg <- end
  }
  mcoef
}


ExpandAllCoefBasis<-function(coef, ncol, value=1, ...){
# coef are the vectorised coef of basis of ncol variables
# coef contains the nbase-1 spline coef for each vars (without intercept)
# add a line of "value"  to the matrix coef
  stopifnot(is.integer(ncol))
  if(length(value)!=ncol){
    value <- rep(value, length.out=ncol)
  }
  mcoef <- matrix(coef, ncol=ncol)
  mcoef <- rbind(value, mcoef)
  mcoef
}


oldExpandAllCoefBasis<-function(coef, ncol, value=1, ...){
# coef are the vectorised coef of basis of ncol variables
# coef contains the nbase-1 spline coef for each vars (without intercept)
# add a line of "value"  to the matrix coef
  stopifnot(is.integer(ncol))
  if(!is.vector(value)){
    value <- rep(value, ncol)
  }
  else {
    stopifnot(length(value)==ncol)
  }
  mcoef <- matrix(coef, ncol=ncol)
  mcoef <- rbind(value, mcoef)
  mcoef
}


