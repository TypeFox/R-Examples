##
## Most Diversified Portfolio
##
PMD <- function(Returns, percentage = TRUE, ...){
  if(is.null(dim(Returns))){
    stop("Argument for 'Returns' must be rectangular.\n")
  }
  call <- match.call()  
  V <- cov(Returns, ...)
  C <- cov2cor(V)
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
  opt <- solve.QP(Dmat = 2 * C, dvec = Dvec, Amat = Amat, bvec = Bvec, meq = meq)
  ## Recovering weights for assets
  w <- opt$solution / sqrt(diag(V))  
  names(w) <- colnames(Returns)
  wnorm <- w / sum(w)
  if(percentage) wnorm <- wnorm * 100
  obj <- new("PortSol", weights = wnorm, opt = opt, type = "Most Diversifified", call = call)
  return(obj)
}
