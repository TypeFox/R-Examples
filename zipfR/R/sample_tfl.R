sample.tfl <- function (obj, N, force.list=FALSE) {
  if (! inherits(obj, "tfl")) stop("first argument must be object of class 'tfl'")
  if (attr(obj, "incomplete")) stop("incomplete type frequency lists are not supported")
  if (! (is.numeric(N) && all(N >= 0))) stop("N must be vector of non-negative integers")
  idx <- order(N)
  if (any(diff(idx) != 1)) {
    warning("N must be increasing (elements have been reordered!)")
    N <- N[idx]
  }
  if (any(N > N(obj))) stop("can't sample ", max(N), " tokens from 'tfl' object with N=", N(obj))
  
  f <- obj$f
  if (! is.integer(f)) {            # convert frequency vector to integer mode
    if (any(f != floor(f))) stop("type frequencies must be integer values in 'obj'")
    f <- as.integer(f)
  }

  result <- list()
  if (length(N) > 0) {
    delta.N <- diff(c(0, N))        # sizes of the consecutive differential samples
    complement <- f                 # frequency vector of complement
    sample <- integer(length(f))    # frequency vector of sample (incremented iteratively)
    for (i in 1:length(N)) {
      split <- zipfR.tflsplit(complement, delta.N[i])
      sample <- sample + split$sample
      complement <- split$complement
      result[[ format(N[i], scientific=FALSE) ]] <-
        tfl(f=sample, k=obj$k, type=obj$type, delete.zeros=TRUE)
    }
  }
    
  if (length(N) == 1 && !force.list) {
    result[[ 1 ]]
  }
  else {
    result
  }
}

## this is the workhorse function: it splits a type frequency vector into two components:
## a) a random sample of specified size N; b) a "complement" list for the remaining tokens
## (no checks for argument types and validity because this function is only used internally)
zipfR.tflsplit <- function (full, N) {
  N <- as.integer(N)
  N.full <- sum(full)
  n.types <- length(full)
  
  if (N > N.full) {
    stop("internal error (N > N(full))")
  }
  else if (N == N.full) {
    list(sample=full, complement=integer(n.types))
  }
  else if (N == 0) {
    list(sample=integer(n.types), complement=full)
  }
  else {
    remain.sample <- N           # number of tokens that still have to be moved to sample
    remain.total <- N.full       # remaining number of tokens that will be considered
    complement <- full           # will become type frequency vector of complement
    sample <- integer(n.types)   # initialized to 0's
    for (k in 1:n.types) {
      ## at this point, we have remain.total tokens left (type IDs k and higher), of which
      ## we need to include remain.sample tokens in our sample; we consider the complement[k]
      ## instances of type ID k and decide 
      n.select <- rhyper(1, remain.sample, remain.total-remain.sample, complement[k])
      remain.total <- remain.total - complement[k]
      ## move n.select instances of type k from complement to sample
      remain.sample <- remain.sample - n.select
      sample[k] <- n.select
      complement[k] <- complement[k] - n.select
    }
    if (remain.sample != 0) stop("internal error (", remain.sample, " missing tokens in sample)")
    list(sample=sample, complement=complement)
  }
}
