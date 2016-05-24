## minimalist partial least squares
pls <- function(x, y, K=1, scale=TRUE, verb=TRUE){

  if(inherits(x,"simple_triplet_matrix")) x <- stm2dg(x)
  xorig <- x
  if(scale){
    scale = sdev(x)
    x <- t(t(x)/scale) 
  }
  else{ scale = NULL }
 
  phi <- matrix(ncol=K, nrow=ncol(x))
  shift <- rep(0,K)
  z <- matrix(ncol=K, nrow=nrow(x))
  yhat <- NULL

  if(ncol(as.matrix(y)) > 1){stop( "PLS only works for univariate y (vector or single column matrix).")}
  v <- scale(as.numeric(y))
  
  if(verb){ cat("Directions ") }
  for(k in 1:K){
    if(verb){ cat(paste(k,", ",sep="")) }

    ## inverse regression, equiv: t(lm(X~v[,k])$coef)[,2]
    phi[,k] <- as.matrix(corr(x,v[,k]))
    
    ## project the fitted direction
    if(inherits(x, "dgCMatrix")){
      z[,k] <- as.matrix(x%*%phi[,k])
    } else { z[,k] <- x%*%phi[,k] }

    ## ortho-normalize
    if(k<K){ v <- cbind(v, lm(v[,k] ~ z[,k])$resid) }
    if(k==1) zlm <- lm(z[,1] ~ 1)
    else zlm <- lm(z[,k] ~ z[,(1:(k-1))])
    s <- sdev(zlm$resid)
    if(s==0) stop( sprintf("perfect fit in pls at k = %d", k) )
    z[,k] <- resid(zlm)/s
    shift[k] <- coef(zlm)[1]/s
    if(k==1) phi[,1] <- phi[,1]/s
    else{
      phi[,k] <- (phi[,k] - phi[,(1:(k-1)),drop=FALSE]%*%coef(zlm)[-1])/s
      shift[k] <- shift[k] - sum(coef(zlm)[-1]*shift[1:(k-1)])/s }

    ## fitted values
    yhat <- cbind(yhat, lm(as.numeric(y)~z[,1:k])$fitted)
  }  
  if(verb){ cat("done.\n")}
  
  fwdmod = lm(as.numeric(y)~z)
  if(!is.null(scale)){ phi <- phi/scale }

  dimnames(z) <- list(rownames(z), direction=paste("z",1:ncol(z),sep=""))
  dimnames(yhat) <- list(rownames(x), model=paste("pls",1:ncol(z),sep=""))
  dimnames(phi) <- list(colnames(x), factor=paste("v",1:ncol(z) -1 ,sep=""))
  names(shift) <- paste("z",1:ncol(z),sep="")
  
  out <- list(y=y, x=xorig, directions=z, loadings=phi, shift=shift,
              fitted = yhat, fwdmod = fwdmod)
  class(out) <- "pls"
  return(out) }


## S3 plot method 
plot.pls <- function(x, K=NULL, xlab="response", ylab=NULL, ...){
  if(is.null(K)){ K <- 1:ncol(x$fitted) }
  y <- x$y 
  K <- as.vector(K)
  par(mfrow=c(1,length(K)))
  ylb <- ylab
  for(k in K){
    if(is.null(ylab)){ ylb=paste("pls(",k,") fitted values") }
    plot(x$fitted[,k] ~y, xlab=xlab, ylab=ylb, ...)
    legend("topleft", legend=paste("corr =",
                        round(cor(as.numeric(y),x$fitted[,k]),2)), bty="n", cex=1.4)
  }
}
  

## S3 method summary function
summary.pls <- function(object, ...){
  print(object)
  cat("Forward regression summary:\n")
  print(summary(object$fwdmod))
}

## S3 method summary function
print.pls <- function(x, ...){
  cat(paste("\nA pls(", ncol(x$directions), ") object, reduced from ", ncol(x$x), " input variables. \n\n", sep="")) }

 ## S3 method predict function
predict.pls <- function(object, newdata, type="response", ...)
{
  if(is.vector(newdata))
    newdata <- matrix(newdata, nrow=1) 
  if(inherits(newdata, "simple_triplet_matrix")) 
    newdata <- stm2dg(newdata)
  if(inherits(newdata,"dgCmatrix"))
    z <- newdata%*%object$loadings
  else z <- as.matrix(newdata)%*%object$loadings
  
  z <- t(t(z) - object$shift)
  if(type=="response"){
    fitted <- cbind(1,z)%*%object$fwdmod$coef
    return(fitted)
  } else{  return(z) }
}

  

