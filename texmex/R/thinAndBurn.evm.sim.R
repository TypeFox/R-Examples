thinAndBurn <- function (object, burn, thin){
  UseMethod("thinAndBurn")
}

thinAndBurn.evmSim <- function(object, burn, thin){

  if(missing(burn)){
    burn <- object$burn
  } else {
    object$burn <- burn
  }
  if(missing(thin)){
    thin <- object$thin
  } else {
    object$thin <- thin
  }
  if(is.null(object$thin)){
    stop("thin or its reciprocal must be a positive integer, for no thinning use thin=1")
  }
  if(is.null(object$burn)){
    stop("burn must be a non-negative integer, for no burn in use burn=0")
  }

  if (thin < 1) thin <- 1 / thin
  if (thin %% 1 > 10^(-6)) stop("thin, or its reciprocal, should be an integer")
  if (burn > dim(object$chains)[1]) stop("burn-in is longer that the whole chain")

  if (burn > 0){
     object$param <- object$chains[-(1:burn), ] # Remove burn-in
  } else {
     object$param <- object$chains
  }
  wh <- 1:nrow(object$param) %% thin == 0
  object$param <- object$param[wh ,]

  invisible(object)
}

