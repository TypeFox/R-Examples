#-----------------------
# Benjamin Risk 
# 25 February 2013
# Edits to David Matteson's PITdCovICA code.
#-----------------------

#---------------------------------------------------
# Benjamin Risk
# 12 March 2013
# modified stats::density.default to return the distribution and the density; 
# note that only the Gaussian kernel works--additional edits needed for other kernels.
#-----------------------------------------------------
BinDist <- function(x,weights,lo,hi,n){
  .Call('BinDistC',x,weights,lo,hi,n)
}

kcdf<-function(x, bw = "SJ", adjust = 1, kernel ="gaussian", weights=NULL, n = 512, from, to, cut = 3) {
  x = sort(x)
  if (!is.numeric(x)) stop("argument 'x' must be numeric")
  x <- as.vector(x)
  x.na <- is.na(x)
  if (any(x.na)) stop("'x' contains missing values")
  N <- nx <- length(x)
  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
    nx <- length(x)
  }
  if (is.null(weights)) {
    weights <- rep.int(1/nx, nx)
    totMass <- nx/N
  }
  else {
    if (length(weights) != N) 
      stop("'x' and 'weights' have unequal length")
    if (!all(is.finite(weights))) 
      stop("'weights' must all be finite")
    if (any(weights < 0)) 
      stop("'weights' must not be negative")
    wsum <- sum(weights)
    if (any(!x.finite)) {
      weights <- weights[x.finite]
      totMass <- sum(weights)/wsum
    }
    else totMass <- 1
    if (!isTRUE(all.equal(1, wsum))) 
      warning("sum(weights) != 1  -- will not get true density")
  }
  n.user <- n
  n <- max(n, 512)
  if (n > 512) 
    n <- 2^ceiling(log2(n))
  if (is.character(bw)) {
    if (nx < 2) 
      stop("need at least 2 points to select a bandwidth automatically")
    bw <- switch(tolower(bw), nrd0 = bw.nrd0(x), nrd = bw.nrd(x), ucv = bw.ucv(x), bcv = bw.bcv(x), sj = , `sj-ste` = bw.SJ(x,method = "ste"), `sj-dpi` = bw.SJ(x, method = "dpi"),stop("unknown bandwidth rule"))
  }
  if (!is.finite(bw)) 
    stop("non-finite 'bw'")
  bw <- adjust * bw
  if (bw <= 0) 
    stop("'bw' is not positive.")
  if (missing(from)) 
    from <- min(x) - cut * bw
  if (missing(to)) 
    to <- max(x) + cut * bw
  if (!is.finite(from)) 
    stop("non-finite 'from'")
  if (!is.finite(to)) 
    stop("non-finite 'to'")
  lo <- from - 4 * bw #They already have from <- min(x) - cut*bw; why this extra?
  up <- to + 4 * bw
  y <- BinDist(x, weights, lo, up, n)*totMass
  kords <- seq.int(0, 2 * (up - lo), length.out = 2L * n)
  kords[(n + 2):(2 * n)] <- -kords[n:2] #What is this doing???
  ##EDITS: original commented out
  #kords <- switch(kernel, gaussian = dnorm(kords, sd = bw))
  kords.temp <- kords
  kords <- pnorm(-kords.temp, sd = bw)
  kords.den <- dnorm(kords.temp, sd = bw)
  rm(kords.temp)
  
  kords <- fft(fft(y) * Conj(fft(kords)), inverse = TRUE)
  kords.den <- fft(fft(y) * Conj(fft(kords.den)), inverse = TRUE)
  kords <- pmax.int(0, Re(kords)[1L:n]/length(y))
  kords.den <- pmax.int(0, Re(kords.den)[1L:n]/length(y))
  
  xords <- seq.int(lo, up, length.out = n)
  #x <- seq.int(from, to, length.out = n.user)
  #y = approx(xords, kords, x)$y
  #rval <- approxfun(x, y,method = "linear", yleft = 0, yright = 1, f = 0, ties = "ordered")
  rval <- approxfun(xords, kords, method = "linear", yleft = 0, yright = 1, f = 0, ties = "ordered")
  denval <- approxfun(xords, kords.den, method = "linear", yleft = lo, yright = up, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  attr(rval, "call") <- sys.call()
  class(denval) <- c("pdf", "fun", class(rval))
  attr(denval, "call") <- sys.call()
  return(list(Fx = rval, fx = denval))
}

#------------------------
est.PIT = function(S, bw='nrd0',adjust = 1){
  n <- nrow(S)
  d <- ncol(S)
  SH = sh = matrix(0,n,d)
  for(j in 1:d){
    KCDF = kcdf(S[,j], bw=bw, adjust = adjust)
    SH[,j] = KCDF$Fx(S[,j])
    sh[,j] = KCDF$fx(S[,j])
  }
  return(list(Fx = SH, fx = sh)) 
}
