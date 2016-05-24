##
## Minimum Tail dependence Portfolio
##
PMTD <- function(Returns, method = c("EmpTC", "EVT"), k = NULL, percentage = TRUE, ...){
  if (is.null(dim(Returns))) {
    stop("Argument for 'Returns' must be rectangular.\n")
  }
  call <- match.call()
  V <- tdc(x = Returns, method = method, k = k, ...)
  N <- ncol(Returns)
  a1 <- rep(1, N)
  b1 <- 1
  a2 <- diag(N)
  b2 <- rep(0, N)
  Amat <- cbind(a1, a2)
  Bvec <- c(b1, b2)
  meq <- c(1, rep(0, N))
  Dvec <- rep(0, N)
  opt <- solve.QP(Dmat = 2 * V, dvec = Dvec, Amat = Amat, bvec = Bvec, meq = meq)
  w <- opt$solution
  ## re-scaling weights by assets' sd 
  sd <- sqrt(diag(cov(Returns)))
  w <- w / sd
  w <- w / sum(w)
  names(w) <- colnames(Returns)
  if(percentage) w <- w * 100
  obj <- new("PortSol", weights = w, opt = opt, type = "Minimum Tail Dependent", call = call)
  return(obj)
}
