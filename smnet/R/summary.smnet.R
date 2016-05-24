summary.smnet<-function(object, ...){
  adjacency  <- object$internals$adjacency
  n.smooth   <- object$internals$n.smooth
  retPrint   <- object$internals$retPrint
  n.linear   <- object$internals$n.linear
  X.list     <- object$internals$X.list
  beta_hat   <- object$internals$beta_hat
  U          <- object$internals$U
  sigma.sq   <- object$internals$sigma.sq
  XTX.spam   <- object$internals$XTX.spam
  X.spam     <- object$internals$X.spam
  lin.names  <- object$internals$lin.names
  ED         <- object$internals$ED
  fit        <- object$internals$fit
  n          <- length(object$internals$response)
  response   <- object$internals$response
  
  # Summarise the linear part of the model
  linearExists <- n.linear + 1
  if(linearExists > 0){
    getcol    <- function(M) ifelse(ncol(M) == "NULL", 1, ncol(M)) 
    cov.dims  <- lapply(X.list, getcol)    
    inds      <- unlist(cov.dims)
    cum.inds  <- cumsum(inds)
    n.cov     <- length(X.list)
    unit.locations <- which(inds == 1)
    
    # get standard errors 
    left1         <- forwardsolve.spam(U, t(X.spam))
    left2         <- backsolve.spam(U, left1)
    separs        <- sqrt(sigma.sq*rowSums(left2^2)[unit.locations])
    pars          <- beta_hat[unit.locations]
    ret.mat       <- cbind(pars, separs, pars/separs, pt(abs(pars/separs), df = n - ED, lower.tail = FALSE) * 2)
    ret           <- as.data.frame(matrix(prettyNum(ret.mat, digits = 5), nrow = length(pars), byrow = F))
    rownames(ret) <- c("(Intercept)", lin.names) 
    colnames(ret) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    cat("\n-----------------------------------------------\nLinear terms:\n-----------------------------------------------\n")
    print(ret)
  }
  if(n.smooth + !is.null(adjacency) > 0){
    cat("\n\n-----------------------------------------------\nSmooth terms:\n-----------------------------------------------\n")
    print(as.data.frame(retPrint))
    cat("\n")
  }
  invisible(list(linear.terms = ret, smooth.terms = as.data.frame(retPrint)))
}




  