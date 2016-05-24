predict.bkpc <-
function(object,  newdata = NULL, n.burnin = 0, ...){


  
  
  if (inherits(object,  "bkpc.kern")){
    if (inherits(object,  "bkpc.kernelMatrix")){
      if(is.null(newdata))testK <- object$x
      else if (inherits(newdata,  "kernelMatrix")) testK <- newdata
      else stop("error: newdata should be of class 'kernelMatrix'") 
    }
    else{
      if(is.null(newdata))testK <- object$x
      else{ 
        if (inherits(newdata,  "kern")) testK <- newdata
        else stop("error: newdata should be of class 'kern'") 
        }
    }
  }
  
  if (!inherits(object,  "bkpc.kern")){
    if(is.null(newdata))testK <- gaussKern(object$x, object$x, object$theta)$K
    else{
      if (inherits(newdata,  "kern"))stop("error: newdata is of class 'kern'") 
      else if (inherits(newdata,  "kernelMatrix"))stop("error: newdata is of class 'kernelMatrix'") 
      else testK <- gaussKern(object$x, newdata, object$theta)$K
    }
  }
  
  
  n.class <- object$n.class
  beta <- as.matrix(object$beta)
  n.samples <- dim(beta)[1]
  n.kpc <- object$n.kpc
  n <- dim(testK)[1]
  
  if (n.burnin >= n.samples) stop("error: too many burn-in iterations specified")
  if (n.burnin < 0) n.burnin <- 0
  
  
  if (object$rotate == TRUE){
    if (object$intercept == TRUE){
      design <- predict(object$kPCA, testK)[, 1 : (n.kpc - 1),  drop = FALSE] 
      design <- cbind(matrix(1, n, 1), design)  
    }
    else design <- predict(object$kPCA, testK)[, 1 : (n.kpc),  drop = FALSE] 
  }
  else{
      design <- testK    
      if (object$intercept == TRUE) design <- cbind(matrix(1, n, 1), design)  
  }
  
  
  e1 <- predictMultinomSamples(design, beta, n.class = n.class, n.burnin)

 # class(e1) <- "bkpcPrediction"
  return(e1)
}
