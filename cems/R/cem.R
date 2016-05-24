


cem <- function (y, x, knnX = 50, sigmaX=  1/3, iter = 100, nPoints = nrow(y),
    stepX = 0.25, stepBW = 0.1, verbose=1, risk=2, penalty = 0, sigmaAsFactor=T, 
    optimalSigmaX = F, quadratic=F ) 
{
    this.call <- match.call()
    if(is.null(nrow(y))){
      y <- as.matrix(y, ncol=1)
    }
    else{ 
      y <- as.matrix(y)
    }    
    if(is.null(nrow(x))){
      x <- as.matrix(x, ncol=1)
    }
    else{ 
      x <- as.matrix(x)
    }    


    nry <- nrow(y)
    nrx <- nrow(x)
    if(nrx != nry){
      stop("coordinate mapping misspecified: y and x don't have the same number of observations") 
    }
    ncy <- ncol(y)
    ncx <- ncol(x)
    res <- .Call("cem_create", as.double(t(y)), nry, ncy,
        as.double(t(x)), nrx, ncx, as.integer(knnX), as.double(sigmaX),
        as.integer(iter), as.integer(nPoints), as.double(stepX),
        as.double(stepBW), as.integer(verbose), as.integer(risk),
        as.integer(penalty), as.integer(sigmaAsFactor),
        as.integer(optimalSigmaX), as.integer(quadratic) )  
   
    obj <- structure( list( y=y, x=t(as.matrix(res[[1]])), knnX =
          knnX, sigmaX = res[[2]], risk=risk,
          penalty = penalty, quadratic=quadratic ), class="cem") 
    obj

}




#Optimize existing CEM further
cem.optimize <- function(object, iter = 100, nPoints = nrow(object$y),
    stepX=1, stepBW=0.1, verbose=1, optimalSigmaX =  F ){

  y <- object$y
  x <- object$x
  knnX <- object$knnX
  sigmaX <- object$sigmaX
  nry <- nrow(y)
  nrx <- nrow(x)
  ncy <- ncol(y)
  ncx <- ncol(x)
  res <- .Call("cem_optimize", as.double(t(y)), nry, ncy, as.double(t(x)),
       nrx, ncx, as.integer(knnX),  as.double(sigmaX), as.integer(iter),
       as.integer(nPoints), as.double(stepX), as.double(stepBW),
       as.integer(verbose), as.integer(object$risk), as.integer(object$penalty),
       as.integer(optimalSigmaX), as.integer(object$quadratic) ) 
    
  object$x = t(as.matrix(res[[1]]))
  object$sigmaX = res[[2]] 
  object
  
}





#Compute geodesic
cem.geodesic <- function(object, xs, xe, iter = 100, step = 0.01,
    verbose=1, ns=100){

  y <- object$y
  x <- object$x
  knnX <- object$knnX
  sigmaX <- object$sigmaX
  nry <- nrow(y)
  nrx <- nrow(x)
  ncy <- ncol(y)
  ncx <- ncol(x)
  res <- .Call("cem_geodesic", as.double(t(y)), nry, ncy, as.double(t(x)), nrx,
      ncx, as.integer(knnX), as.double(sigmaX), as.integer(iter),
      as.double(step), as.integer(verbose), as.integer(object$quadratic),
      as.double(xs), as.double(xe), as.integer(ns)) 
    
  t(res)
}





predict.cem <- function(object, newdata = object$y, type=c("coordinates",
      "curvature" ), ... ){

  type=match.arg(type)

  if( is.null( nrow(newdata) ) ){
    data <- as.matrix(newdata, ncol=1)
  }
  else{
    data <- as.matrix(newdata)
  }
  y <- object$y
  x <- object$x
  knnX <- object$knnX
  sigmaX <- object$sigmaX
  nry <- as.integer(nrow(y))
  nrx <- as.integer(nrow(x))
  ncy <- as.integer(ncol(y))
  ncx <- as.integer(ncol(x))
  nrd <- as.integer(nrow(data))

  if(type == "coordinates"){

    if(ncol(data)  == ncol(object$y) ){
      res <- .Call("cem_parametrize", as.double(t(data)), nrd,
          as.double(t(y)), nry, ncy, as.double(t(x)), nrx, ncx,
          as.integer(knnX), as.double(sigmaX), as.integer(object$quadratic)
          )  
      res <- t(as.matrix(res))
    }
    else{ 
      res <- .Call("cem_reconstruct", as.double(t(data)), nrd,
          as.double(t(y)), nry, ncy, as.double(t(x)), nrx, ncx,
          as.integer(knnX), as.double(sigmaX), as.integer(object$quadratic))
      tangents = c()
      for(i in 1:ncx){
        tangents[[i]] = t(as.matrix(res[[1+i]]))
      }
      res <- list(y = t(as.matrix(res[[1]])), tangents=tangents)
    }
  }

  else if(type == "curvature"){
     res <- .Call("cem_curvature", as.double(t(data)), nrd,
     as.double(t(y)), nry, ncy, as.double(t(x)), nrx, ncx, as.integer(knnX),
     as.double(sigmaX), as.integer(object$quadratic))
     
       res <- list( principal = t(as.matrix(res[[1]])), 
                   mean =   as.vector(res[[2]]), 
                   gauss =  as.vector(res[[3]]), 
                   detg =   as.vector(res[[3]]), 
                   detB =   as.vector(res[[4]]), 
                   frob =   as.vector(res[[5]]) 
                 )
     
  }

res
}


