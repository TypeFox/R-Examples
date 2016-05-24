# Sensitivity indices based on kernel embeddings of distributions
#
# Sebastien Da Veiga 2014


# Kernel functions
rbf_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(exp(-0.5*d^2))
}

laplace_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(exp(-d))
}

dcov_hsic <- function(x,param,d=NULL){
  nobs <- length(x)
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  return(0.5*(-d^param+matrix(abs(x),nobs,nobs)^param+t(matrix(abs(x),nobs,nobs))^param))
}

raquad_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(1-d^2/(1+d^2))
}

invmultiquad_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return(1/sqrt(d^2+1))
}

linear_hsic <- function(x,...){
  return(x%*%t(x))
}

matern3_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return((1+sqrt(3)*d)*exp(-sqrt(3)*d))
}

matern5_hsic <- function(x,param,d=NULL){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  d <- d/param
  return((1+sqrt(5)*d+5/3*d^2)*exp(-sqrt(5)*d))
}

# Bernoulli polynomials
B1 <- function(t){
  return(t-.5)
}
B2 <- function(t){
  return(t^2-t+1/6)
}
B4 <- function(t){
  return(t^4-2*t^3+t^2-1/30)
}

ssanova1_hsic <- function(x,d=NULL,...){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  return(1+B1(x)%*%t(B1(x))+0.5*B2(d))
}

ssanova2_hsic <- function(x,d=NULL,...){
  if (is.null(d)){
    d <- as.matrix(dist(x))
  }
  return(1+B1(x)%*%t(B1(x))+B2(x)%*%t(B2(x))/4-B4(d)/24)
}


HSIC <- function(X, Y, kernelX, paramX, kernelY, paramY){
  
  nobs <- nrow(X)
  
  # Centering matrix
  H <- diag(nobs) - 1/nobs*matrix(1,nobs,nobs)
  
  KX <- do.call(get(paste(kernelX,"_hsic",sep="")), list(x=X,param=paramX))
  KXH <- H%*%KX%*%H
  
  KY <- do.call(get(paste(kernelY,"_hsic",sep="")), list(x=Y,param=paramY))
  KYH <- H%*%KY%*%H
  
  return(list(estimate=mean(KXH*KYH),meanX=mean(KXH*KXH),meanY=mean(KYH*KYH)))
}


sensiHSIC <- function(model = NULL, X, kernelX = "rbf", paramX = NA, 
                      kernelY = "rbf", paramY = NA, nboot = 0, conf = 0.95, ...) {
  
  
  if (is.data.frame(X)){
    X <- as.matrix(unname(X))
  }
  
  else if(!is.matrix(X)){
    stop("The sample X must be a matrix or a data frame")
  }
  
  p <- ncol(X)
  nkx <- length(kernelX)
  if (!(nkx == 1 | nkx == p)){
    stop("KernelX must be of length 1 or p (number of input variables)")
  }
  if (!missing(paramX)){
    npx <- length(paramX)
    if (!(npx == 1 | npx == p)){
      stop("paramX must be of length 1 or p (number of input variables)")
    }
  }
  
  x <- list(model = model, X = X, kernelX = kernelX, paramX = paramX, 
            kernelY = kernelY, paramY = paramY, nboot = nboot,
            conf = conf, call = match.call()) 
  class(x) <- "sensiHSIC"
  
  #calcul of the response for explicit model
  if (! is.null(x$model)) {
    response(x, ...)
    x=tell(x, ...)
  }  
  return(x)
}


estim.sensiHSIC <- function(data, i=1:nrow(data), kernelX, paramX, 
                            kernelY, paramY) {
  
  ptot <- ncol(data)
  p <- ptot - 1
  X <- data[i,1:p]
  Y <- data[i,ptot]
  S = matrix(0,nrow=p,ncol=1)
  
  # HSIC indices
  for (i in 1:p){
    Xtemp <- as.matrix(X[,i])
    Ytemp <- as.matrix(Y)
    res <- HSIC(Xtemp, Ytemp, kernelX[i], paramX[i], kernelY, paramY)
    S[i] <- res$estimate/sqrt(res$meanX*res$meanY)
  }
  return(S)
}


tell.sensiHSIC <- function(x, y = NULL, ...) {

  id <- deparse(substitute(x))
  if (! is.null(y)) {
    x$y <- y
  } 
  else if (is.null(x$y)) {
    stop("y not found")
  }
  
  n <- nrow(x$X)
  p <- ncol(x$X)
  
  nkx <- length(x$kernelX)
  if (nkx==1) x$kernelX <- rep(x$kernelX,p)
  
  if (is.na(x$paramX[1])){
    x$paramX <- matrix(nrow=1,ncol=p)
    for (i in 1:p){
      if (x$kernelX[i]=="dcov"){
        x$paramX[i] = 1
      }
      else{
        x$paramX[i] = sd(x$X[,i])
      }
    }
  }
  else{
    if (length(x$paramX)==1) x$paramX <- rep(x$paramX,p)
  }
  
  
  if (is.na(x$paramY)){
    if (x$kernelY=="dcov"){
      x$paramY = 1
    }
    else{
      x$paramY <- sd(x$y)
    }
  }
  
  data <-cbind(x$X,x$y)
  
  if (x$nboot == 0) {
    res <- estim.sensiHSIC(data, 1:n, x$kernelX, x$paramX, x$kernelY, x$paramY)
    x$S <- data.frame(res)
    colnames(x$S) <- "original"
  } 
  else {
    S.boot <- boot(data, estim.sensiHSIC, 
                   kernelX = x$kernelX, paramX = x$paramX, 
                   kernelY = x$kernelY, paramY = x$paramY, R = x$nboot)
    x$S <- bootstats(S.boot, x$conf, "basic")
  }
  
  
  rownames <- c()
  for (i in 1:p) {
    rownames <- c(rownames,paste("S",i,sep=""))
  }
  rownames(x$S) <- rownames
  
  assign(id, x, parent.frame())
  return(x)
}


print.sensiHSIC<- function(x, ...) {
  
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    if (! is.null(x$S)) {
      cat("\n\n\nHSIC indices\n")
      print(x$S)
    }
    else{
      cat("(empty)\n")
    }
  }
}


plot.sensiHSIC <- function(x, ylim = c(0, 1), ...) {
  
  if (! is.null(x$y)) {
    nodeplot(x$S, ylim = ylim)
    legend(x = "topright", legend = "HSIC indices")
  }
}
