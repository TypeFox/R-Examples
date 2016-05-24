"stepsel" <-
function(path, y, type = c("MRC", "AIC", "BIC"), ...)
{
  type <- match.arg(type)
  switch(type,
    MRC = stepsel.MRC(path, y, ...),
    AIC = stepsel.AIC(path, ...),
    BIC = stepsel.BIC(path, ...),
    stop("type ", type, " not known!")
  )
}

"stepsel.MRC" <-
function(path, y, q, alpha = 0.05, r = ceiling(50 / min(alpha, 1 - alpha)), lengths = if(attr(path$cand, "family") == "gaussKern") 2^(floor(log2(length(y))):ceiling(log2(length(attr(path$cand, "param")$kern)))) else 2^(floor(log2(length(y))):0), penalty = c("none", "log", "sqrt"), name = if(attr(path$cand, "family") == "gaussKern") ".MRC.ktable" else ".MRC.table", pos = .GlobalEnv)
{
  family <- attr(path$cand, "family")
  stopifnot(family %in% c("gauss", "gaussKern"))
  if(is.character(penalty)) {
    penalty <- match.arg(penalty)
  }
  
  # data and kernel lengths
  n <- length(y)
  kl <- if(family == "gaussKern") length(attr(path$cand, "param")$kern) else 1
  
  # precompute
  chiFFT <- chi.FFT(n, lengths)
  # valid intervals
  K <- sapply(n - lengths + 1, function(k) rep(c(T,F), c(k, n - k)))
  
  # intialise
  MRCs <- rep(NA, length(path))
  m <- 0
  val <- Inf
  
  # compute quantile
  if(missing(q)) {
    if(is.null(r)) stop("q or r need to be specified!")
    q <- if(family == "gaussKern") {
      kMRC.quant(1 - alpha, n, r, attr(path$cand, "param")$kern, lengths, name, pos)
    } else {
      MRC.quant(1 - alpha, n, r, lengths, penalty, name, pos)
    }
  }

  # loop until hypothesis is not rejected
  while(val > q && m < length(path)) {
    m <- m + 1
    fit <- path[[m]]
    res <- resid(path[[m]], y)
    res <- res / sd(res)
    if(family == "gaussKern") res[neighbours(fit$rightEnd, r = kl)] <- 0
    MRCs[m] <- val <- MRC.FFT(as.vector(mvfft(as.matrix(res))), chiFFT, K = K, lengths = lengths, penalty = penalty)["max"]
  }
  if(val > q) warning("Multiresolution criterion never fulfilled!")
  structure(sum(!is.na(MRCs)), crit = MRCs)
}

"stepsel.AIC" <-
function(path, ...)
{
  crit = AIC(path)
  structure(which.min(crit), crit = crit)
}

"stepsel.BIC" <-
function(path, ...)
{
  crit = BIC(path)
  structure(which.min(crit), crit = crit)
}
