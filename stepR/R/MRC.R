## MRC criterion
# compute test intervals with given lengths, dyadic by default (longest first)
chi <- function(n, lengths = 2^(floor(log2(n)):0), norm = sqrt(lengths)) {
  # normalised test intervals (not centered around 0!)
  sapply(lengths, function(k) rep(1:0, c(k, n - k))) / rep(norm, each = n)
}

# compute FFT of test intervals with given lengths, dyadic by default (longest first)
chi.FFT <- function(n, lengths = 2^(floor(log2(n)):0), norm = sqrt(lengths)) {
  # normalised test intervals (not centered around 0!)
  chi <- sapply(lengths, function(k) rep(1:0, c(k, n - k))) / rep(norm, each = n)
  mvfft(chi)
}

# compute MR coefficients using precomputed FFTs, epsFFT must be a vector, i.e. represent a single data set, K is the set of valid intervals
MRCoeff.FFT <- function(epsFFT, testFFT, K = matrix(TRUE, nrow(testFFT), ncol(testFFT)), signed = FALSE) {
  val <- Re(mvfft(epsFFT * Conj(testFFT), inverse = TRUE) / length(epsFFT))
  if(!signed) val <- abs(val)

  # restrict to valid set K
  val[!K] <- NA
  val
}

# compute MR coefficients
MRCoeff <- function(x, lengths = 2^(floor(log2(length(x))):0), norm = sqrt(lengths), signed = FALSE) {
  n <- length(x)
  n2 <- 2^ceiling(log2(n)) # bring to dyadic length
  K <- sapply(n - lengths + 1, function(k) rep(c(TRUE, FALSE), c(k, n2 - k)))
  MRCoeff.FFT(fft(c(x, rep(0, n2 - n ))), chi.FFT(n2, lengths, norm = norm), K = K, signed = signed)[1:n,]
}

# compute MRC using precomputed FFTs, epsFFT must be a vector, i.e. represent a single data set, K is the set of valid intervals
MRC.FFT <- function(epsFFT, testFFT, K = matrix(TRUE, nrow(testFFT), ncol(testFFT)), lengths, penalty = c("none", "log", "sqrt")) {
  if(is.character(penalty)) {
    penalty <- match.arg(penalty)
  }
  n <- length(epsFFT)
  n2 <- nrow(testFFT)
  
  # compute MR coefficients
  val <- if(is.character(penalty)) {
    switch(penalty,
      none =MRCoeff.FFT(epsFFT, testFFT = testFFT, K = K)^2 / 2,
      log = MRCoeff.FFT(epsFFT, testFFT = testFFT, K = K)^2 / 2 + rep( log(lengths / n), each = n2),
      sqrt = abs(MRCoeff.FFT(epsFFT, testFFT = testFFT, K = K)) - sqrt(2 * ( 1 - rep( log(lengths / n), each = n2) )),
      stop("unknown penalty")
    )
  } else {
    penalty(MRCoeff.FFT(epsFFT, testFFT = testFFT, K = K), lengths)
  }
  
  # find maximal value
  index <- which.max(val)[1]
  
  # return position of maximum, corresponding interval length and max. value
  c(ind.x = ( index - 1) %% n + 1, ind.d = ( index - 1 ) %/% n + 1, max = val[index])
}

# compute MRC of signal x
MRC <- function(x, lengths = 2^(floor(log2(length(x))):0), norm = sqrt(lengths), penalty = c("none", "log", "sqrt")) {
  n <- length(x)
  n2 <- 2^ceiling(log2(n)) # bring to dyadic length
  K <- sapply(n - lengths + 1, function(k) rep(c(TRUE, FALSE), c(k, n2 - k)))
  MRC.FFT(fft(c(x, rep(0, n2 - n ))), chi.FFT(n2, lengths, norm = norm), K = K, lengths = lengths, penalty = penalty)
}

# simulate MRC values with given number of repititions r and lengths, dyadic by default (longest first)
MRC.simul <- function(n, r, lengths = 2^(floor(log2(n)):0), penalty = c("none", "log", "sqrt")) {
  if(is.character(penalty)) {
    penalty <- match.arg(penalty)
  }

  # precompute
  n2 <- 2^ceiling(log2(n)) # bring to dyadic length
  chiFFT <- chi.FFT(n2, lengths)
  K <- sapply(n - lengths + 1, function(k) rep(c(TRUE, FALSE), c(k, n2 - k)))
  
  # simulate
  simul <- sapply(1:r, function(i) {
    
    eps <- c(rnorm(n), rep(0, n2 - n))
    epsFFT <- as.vector(mvfft(as.matrix(eps)))
    
    # compute MRCs
    if(is.character(penalty)) {
      switch(penalty,
	none = max(MRCoeff.FFT(epsFFT, testFFT = chiFFT, K = K)^2 / 2, na.rm = TRUE),
	log = max(MRCoeff.FFT(epsFFT, testFFT = chiFFT, K = K)^2 / 2 + rep( log(lengths / n), each = n2), na.rm = TRUE),
	sqrt = max(abs(MRCoeff.FFT(epsFFT, testFFT = chiFFT, K = K)) - sqrt(2 * ( 1 - rep( log(lengths / n), each = n2) )), na.rm = TRUE),
	stop("unknown penalty")
      )
    } else {
      max(penalty(MRCoeff.FFT(epsFFT, testFFT = chiFFT, K = K), lengths), na.rm = TRUE)
    }
  })
    
  sort(simul)
}

# compute MRC p-value(s) for sequence of length n evaluated at lengths based on at least r repetitions, first checking in a table whether this has been simulated before
MRC.pvalue <- function(q, n, r, lengths = 2^(floor(log2(n)):0), penalty = c("none", "log", "sqrt"), name = ".MRC.table", pos = .GlobalEnv, inherits = TRUE) {
  if(is.character(penalty)) {
    penalty <- match.arg(penalty)
    pen <- penalty
  } else if(is.function(penalty)) {
    pen <- deparse(penalty)
  } else {
    stop("penalty must be a character or function")
  }
  index <- 0
  if(exists(name, where = pos, inherits = inherits)) {
    tab <- get(name, pos = pos, inherits = inherits)
    # check whether this has been simulated yet
    check <- sapply(1:length(tab$n), function(i)
      if(tab$n[i] == n & length(tab$lengths[[i]]) == length(lengths)) all(tab$lengths[[i]] == lengths) & all(tab$penalty[[i]] == pen) else F)
    if(any(check))
      index <- which(check)
    else {
      # create entries
      index <- length(tab$n) + 1
      tab$n[index] <- n
      tab$r[index] <- 0
      tab$lengths[[index]] <- lengths
      tab$penalty[[index]] <- penalty
      tab$simul[[index]] <- NA
    }
  } else {
    # create table
    index <- 1
    tab <- list(n = n, r = 0, lengths = list(lengths), penalty = list(pen), simul = list(NA))
  }
  # check whether precomputed simulations suffice
  if(tab$r[index] >= r)
    simul <- tab$simul[[index]]
  else {
    # simulate
    simul <- MRC.simul(n, r, lengths, penalty)
    tab$simul[[index]] <- simul
    tab$r[index] <- r
    # save table
    assign(name, tab, pos = pos, inherits = inherits)
  }
  # return p-value(s) from upper tail
  pmin(1 - ecdf(simul)(q) + 1/length(simul), 1)
}

# compute MRC quantile for sequence of length n evaluated at lengths based on at least r repetitions, first checking in a table whether this has been simulated before
MRC.quant <- function(p, n, r, lengths = 2^(floor(log2(n)):0), penalty = c("none", "log", "sqrt"), name = ".MRC.table", pos = .GlobalEnv, inherits = TRUE, ...) {
  if(is.character(penalty)) {
    penalty <- match.arg(penalty)
    pen <- penalty
  } else if(is.function(penalty)) {
    pen <- deparse(penalty)
  } else {
    stop("penalty must be a character or function")
  }
  index <- 0
  if(exists(name, where = pos, inherits = inherits)) {
    tab <- get(name, pos = pos, inherits = inherits)
    # check whether this has been simulated yet
    check <- sapply(1:length(tab$n), function(i)
      if(tab$n[i] == n & length(tab$lengths[[i]]) == length(lengths)) all(tab$lengths[[i]] == lengths) & all(tab$penalty[[i]] == pen) else F)
    if(any(check))
      index <- which(check)
    else {
      # create entries
      index <- length(tab$n) + 1
      tab$n[index] <- n
      tab$r[index] <- 0
      tab$lengths[[index]] <- lengths
      tab$penalty[[index]] <- penalty
      tab$simul[[index]] <- NA
    }
  } else {
    # create table
    index <- 1
    tab <- list(n = n, r = 0, lengths = list(lengths), penalty = list(pen), simul = list(NA))
  }
  # check whether precomputed simulations suffice
  if(tab$r[index] >= r)
    simul <- tab$simul[[index]]
  else {
    # simulate
    simul <- MRC.simul(n, r, lengths, penalty)
    tab$simul[[index]] <- simul
    tab$r[index] <- r
    # save table
    assign(name, tab, pos = pos, inherits = inherits)
  }
  # return quantile
  quantile(x = simul, probs = p, ...)
}

# simulate MRC values with given number of repititions r and lengths and kernel
kMRC.simul <- function(n, r, kern, lengths = 2^(floor(log2(n)):ceiling(log2(length(kern))))) {

  # precompute
  n2 <- 2^ceiling(log2(n + 2 * length(kern) - 2)) # bring to dyadic length, needs to be padded because of convolution
  chiFFT <- chi.FFT(n2, lengths)
  kernFFT <- Conj(as.vector(mvfft(as.matrix(c(rev(kern), rep(0, n2 - length(kern)))))))
  
  # valid intervals
  K <- sapply(n, function(k) rep(c(TRUE,FALSE), c(k, n2 - k)))
#   K <- rbind(matrix(F, length(kern) - 1, ncol(K)), K, matrix(F, length(kern) - 1, ncol(K)))

  # kernel reduces standard deviation, so we need to start with this
  ksd <- 1 / sqrt(sum(kern^2))
  
  # simulate
  simul <- sapply(1:r, function(i) {
    
    eps <- rnorm(n2, sd = ksd)
    epsFFT <- as.vector(mvfft(as.matrix(eps))) * kernFFT
    
    # compute MRCs
    max(MRCoeff.FFT(epsFFT, testFFT = chiFFT, K = K)^2 / 2, na.rm = TRUE)
  })
    
  sort(simul)
}
# compute MRC p-value(s) for sequence of length n evaluated at lengths based on at least r repetitions, first checking in a table whether this has been simulated before
# lengths start at kernel k's length
kMRC.pvalue <- function(q, n, r, kern, lengths = 2^(floor(log2(n)):ceiling(log2(length(kern)))), name = ".MRC.ktable", pos = .GlobalEnv, inherits = TRUE) {
  index <- 0
  if(exists(name, where = pos, inherits = inherits)) {
    tab <- get(name, pos = pos, inherits = inherits)
    # check whether this has been simulated yet
    check <- sapply(1:length(tab$n), function(i) if(tab$n[i] == n & length(tab$lengths[[i]]) == length(lengths) & length(tab$kern[[i]]) == length(kern))
      all(tab$lengths[[i]] == lengths) && all(tab$kern[[i]] == kern) else F)
    if(any(check))
      index <- which(check)
    else {
      # create entries
      index <- length(tab$n) + 1
      tab$n[index] <- n
      tab$r[index] <- 0
      tab$lengths[[index]] <- lengths
      tab$kern[[index]] <- kern
      tab$simul[[index]] <- NA
    }
  } else {
    # create table
    index <- 1
    tab <- list(n = n, r = 0, lengths = list(lengths), kern = list(kern), simul = list(NA))
  }
  # check whether precomputed simulations suffice
  if(tab$r[index] >= r)
    simul <- tab$simul[[index]]
  else {
    # simulate
    simul <- kMRC.simul(n, r, kern, lengths)
    tab$simul[[index]] <- simul
    tab$r[index] <- r
    # save table
    assign(name, tab, pos = pos, inherits = inherits)
  }
  # return p-value(s) from upper tail
  pmin(1 - ecdf(simul)(q) + 1/length(simul), 1)
}
# compute MRC quantile(s) for sequence of length n evaluated at lengths based on at least r repetitions, first checking in a table whether this has been simulated before
# lengths start at kernel k's length
kMRC.quant <- function(p, n, r, kern, lengths = 2^(floor(log2(n)):ceiling(log2(length(kern)))), name = ".MRC.ktable", pos = .GlobalEnv, inherits = TRUE, ...) {
  index <- 0
  if(exists(name, where = pos, inherits = inherits)) {
    tab <- get(name, pos = pos, inherits = inherits)
    # check whether this has been simulated yet
    check <- sapply(1:length(tab$n), function(i) if(tab$n[i] == n & length(tab$lengths[[i]]) == length(lengths) & length(tab$kern[[i]]) == length(kern))
      all(tab$lengths[[i]] == lengths) && all(tab$kern[[i]] == kern) else F)
    if(any(check))
      index <- which(check)
    else {
      # create entries
      index <- length(tab$n) + 1
      tab$n[index] <- n
      tab$r[index] <- 0
      tab$lengths[[index]] <- lengths
      tab$kern[[index]] <- kern
      tab$simul[[index]] <- NA
    }
  } else {
    # create table
    index <- 1
    tab <- list(n = n, r = 0, lengths = list(lengths), kern = list(kern), simul = list(NA))
  }
  # check whether precomputed simulations suffice
  if(tab$r[index] >= r)
    simul <- tab$simul[[index]]
  else {
    # simulate
    simul <- kMRC.simul(n, r, kern, lengths)
    tab$simul[[index]] <- simul
    tab$r[index] <- r
    # save table
    assign(name, tab, pos = pos, inherits = inherits)
  }
  # return quantile
  quantile(x = simul, probs = p, ...)
}
