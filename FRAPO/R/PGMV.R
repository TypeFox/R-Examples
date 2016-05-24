##
## Global Minimum Variance Portfolio
##
PGMV <- function(Returns, percentage = TRUE, ...){
  if(is.null(dim(Returns))){
    stop("Argument for 'Returns' must be rectangular.\n")
  }
  call <- match.call()  
  V <- cov(Returns, ...)
  N <- ncol(Returns)
  ## QP
  ## Budget constraint (equality)
  a1 <- rep(1, N)
  b1 <- 1
  ## Nonnegativity constraint (inequality)
  a2 <- diag(N)
  b2 <- rep(0, N)
  ## combining restrictions
  Amat <- cbind(a1, a2)
  Bvec <- c(b1, b2)
  meq <- c(1, rep(0, N))
  Dvec <- rep(0, N)
  ## Call to solver
  opt <- solve.QP(Dmat = 2 * V, dvec = Dvec, Amat = Amat, bvec = Bvec, meq = meq)
  ## Recovering weights for assets
  w <- opt$solution
  names(w) <- colnames(Returns)
  if(percentage) w <- w * 100
  obj <- new("PortSol", weights = w, opt = opt, type = "Global Minimum Variance", call = call)
  return(obj)
}
