##
## Equal Risk Contribution Portfolio
##
PERC2 <- function (Sigma, par = NULL, percentage = TRUE, ...){
  if(!isSymmetric(Sigma)){
    stop("Matrix provided for Sigma is not symmetric.\n")
  }
  N <- ncol(Sigma)
  if(is.null(par)){
    par <- rep(1/N, N)
  } else {
    if(length(par) != N){
      stop("Length of 'par' not comformable with dimension of 'Sigma'.\n")
    }
  }
  call <- match.call()
  ## objective
  f <- function(pars, Sigma){
    n <- ncol(Sigma)
    rowprod <- c(Sigma %*% pars)
    t1 <- n * sum(pars^2 * (rowprod)^2)
    t2 <- sum(c(outer(pars * rowprod, pars * rowprod)))
    return(t1 - t2)
  }
  opt <- solnp(pars = par, fun = f, LB = rep(0, N), UB = rep(1, N), Sigma = Sigma, ...) 
  w <- opt$par
  w <- w / sum(w)
  if(percentage) w <- w * 100
  if(is.null(dimnames(Sigma))){
    names(w) <- paste("Asset", 1:N, sep = "")
    } else {
      names(w) <- colnames(Sigma)
    }
  obj <- new("PortSol", weights = w, opt = opt,
             type = "Equal Risk Contribution", call = call)
  return(obj)
}
