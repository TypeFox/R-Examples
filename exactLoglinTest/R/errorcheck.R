#x is the design matrix
errorcheck <- function(y, x, stat, dens, nosim, method, savechain, tdf, maxiter, p, batchsize){
  if ((!is.double(y)) & (!is.integer(y)))
    stop("y must be double or integer valued")
  else {
    if (any(y < 0))
      stop("y must be positive")
    if (any(round(y) != y))
      stop("round(y) must equal y")
  }

  if (!is.matrix(x))
    stop("x must a matrix")
  else if ((!is.double(x)) & (!is.integer(x)))
    stop("x must be double or integer valued")
  else if (qr(x)$rank > dim(x)[2])
    stop("Rank(x) <= number of rows")
  else if (dim(x)[2] < 2)
    stop("Model must have more than 2 parameters")
  else if (dim(x)[1] < dim(x)[2])
    stop("Model has more parameters than datapoints")
  else if (qr(x)$rank != dim(x)[2])
    stop("Design matrix must be of full rank")
  
  if (!is.function(stat)) stop("stat must be a function")
  
  if (!is.function(dens)) stop("dens must be a function")

  if ((!is.double(nosim)) & (!is.integer(nosim))) stop("nosim must be double or integer valued")
  else if (nosim <= 0) stop("nosim < 0 not allowed")

  if (method != "bab" & method != "cab") stop("method must be either cab or bab")

  if ((!is.double(tdf)) & (!is.integer(tdf))) stop("tdf must be double or integer valued")
  else if (tdf <= 0) stop("tdf < 0 not allowed")

  if (!is.null(maxiter) & !is.double(maxiter) & !is.integer(maxiter))
    stop("maxiter must be null, double or integer valued")
  else if (method == "bab" & is.null(maxiter)) stop("maxiter must be specified if method = bab")
  else if (maxiter < nosim) stop("maxiter >= nosim")
   
  if (!is.null(p) & !is.double(p) & !is.integer(p))
    stop("p must be null, double or integer valued")
  else if (method == "cab" & is.null(p)) stop("p must be specified if method = cab")
  else if (!is.null(p))
    if (p < 0 | p > 1) stop("p must be between 0 and 1")

  if (!is.null(batchsize) & !is.double(batchsize) & !is.integer(batchsize))
    stop("batchsize must be null, double or integer valued")
  else if (method == "cab"){ 
    if (is.null(batchsize))
      stop("batchsize must be specified if method = bab") 
    else if (batchsize <= 0) stop("batchsize must be > 0")
    else if (nosim / batchsize < 1) stop("batches are too large")
  }
}
