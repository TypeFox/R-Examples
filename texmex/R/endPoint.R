endPoint <- function(y,verbose=TRUE,.unique=TRUE,...){
  UseMethod("endPoint", y)
}

endPoint.evmOpt <- function(y, verbose=TRUE,.unique=TRUE,...){

  if(.unique) Unique <- unique else Unique <- identity

  p <- texmexMakeParams(coef(y), y$data$D)
  endpoint <- y$family$endpoint

  negShape <- p[, ncol(p)] < 0

  if(any(negShape)){
    UpperEndPoint <- endpoint(p, y)
    UpperEndPoint[!negShape] <- Inf
    if(verbose){
      o <- Unique(cbind(y$data$D[['xi']], p))
      print(signif(o,...))
    } else {
      invisible(Unique(UpperEndPoint))
    }
  } else {
      Unique(rep(Inf,length(negShape)))
  }
}

endPoint.evmBoot <- endPoint.evmSim <- function(y,verbose=TRUE,.unique=TRUE,...){
  endPoint(y$map,verbose=verbose,.unique=.unique,...)
}

