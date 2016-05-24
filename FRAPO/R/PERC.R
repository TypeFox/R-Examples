##
## Equal Risk Contribution Portfolio
##
PERC <- function (Sigma, par = NULL, percentage = TRUE, ...){
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
  f <- function(x, Sigma){
    pr <- sqrt(t(x) %*% Sigma %*% x)
    mrc <- c(x * Sigma %*% x) / pr
    val <- sd(mrc)
    return(val)
  }
  opt <- nlminb(start = par, objective = f, Sigma = Sigma,
                lower = 0, upper = 1, ...)
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
