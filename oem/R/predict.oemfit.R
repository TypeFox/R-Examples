predict.oemfit <- function(object, newx, s = NULL,
                           type = c("response",
                             "coefficients",
                             "nonzero"), ...) {
  type <- match.arg(type)
  nbeta <- object$beta
  if(!is.null(s)){
    lambda <- object$lambda
    lamlist <- lambda.interp(object$lambda,s)
    nbeta <- nbeta[,lamlist$left,drop=FALSE]*lamlist$frac +nbeta[,lamlist$right,drop=FALSE]*(1-lamlist$frac)
  }
  if (type == "coefficients") return(nbeta)
  if (type == "nonzero") {
    newbeta <- abs(as.matrix(object$beta)) > 0
    index <- 1:(dim(newbeta)[1])
    nzel <- function(x, index) if(any(x)) index[x] else NULL
    betaList <- apply(newbeta, 2, nzel, index)
    return(betaList)
  }

  newx <- as.matrix(newx)
  # add constant column if needed                                        
  if (ncol(newx) < nrow(nbeta))
    newx <- cbind(rep(1, nrow(newx)), newx)
  nfit <- as.matrix(newx %*% nbeta)
}
