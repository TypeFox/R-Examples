## CF1

cf1 <- function(size, alpha, rate, class="CsparseMatrix") {
  if (missing(size)) {
    if (missing(alpha) || missing(rate)) {
      stop("alpha and rate are needed.")
    } else {
      size <- length(alpha)
    }
  } else {
    if (!missing(alpha) && !missing(rate)) {
      warning("size is ignored.")
      size <- length(alpha)
    } else {
      if (!missing(alpha) || !missing(rate)) {
        warning("alpha and rate are ignored.")
      }
      alpha <- rep(1.0/size, size)
      rate <- rep(1.0, size)      
    }
  }

  xi <- numeric(size)
  xi[size] <- rate[size]
  if (size >= 2) {
    i <- c(1:size, 1:(size-1))
    j <- c(1:size, 2:size)
    x <- c(-rate, rate[1:(size-1)])
    Q <- sparseMatrix(i=i, j=j, x=x)
  } else {
    Q <- matrix(-rate[1],1,1)
  }
  if (!is(Q, class)) {
    Q <- as(Q, class)
  }
  zero <- 1.0e-8
  df <- sum(abs(alpha) > zero) - 1 + sum(abs(rate) > zero)
  new("cf1", size=size, alpha=alpha, Q=Q, xi=xi, rate=rate, df=df)
}

cf1.sample <- function(n, ph) {
  res <- rexp(n, ph@rate[ph@size])
  if (ph@size == 1) {
    return(res)
  }
  x <- cumsum(rmultinom(n=1, size=n, prob=ph@alpha))
  for (l in (ph@size-1):1) {
    y <- x[l]
    if (y == 0) break
    res[1:y] <- res[1:y] + rexp(y, ph@rate[l])
  }
  sample(res)
}

cf1.param <- function(size, data,
  diff.init = c(1, 4, 16, 64, 256, 1024),
  scale.init = c(0.5, 1.0, 2.0), maxiter.init = 5, verbose, class) {
  if (verbose$cf1init) cat("Initializing CF1 ...\n")
  m <- mapfit.mean(data)
  maxllf <- -Inf
  for (x in scale.init) {
    for (s in diff.init) {
      ph <- cf1.param.linear(size, m * x, s, class)
      phres <- try(emfit(model=ph, data=data,
        initialize = FALSE, control = list(maxiter=maxiter.init)), silent = TRUE)
      if (class(phres) != "try-error") {
        if (is.finite(phres$llf)) {
          if (maxllf < phres$llf) {
            maxllf <- phres$llf
            maxph <- ph
            if (verbose$cf1init) cat("o")
          }
          else {
            if (verbose$cf1init) cat("x")
          }
        }
        else {
          if (verbose$cf1init) cat("-")
        }
      }
      else {
        if (verbose$cf1init) cat("-")
      }
    }
    if (verbose$cf1init) cat("\n")
    for (s in diff.init) {
      ph <- cf1.param.power(size, m * x, s, class)
      phres <- try(emfit(model=ph, data=data,
        initialize = FALSE, control = list(maxiter=maxiter.init)), silent = TRUE)
      if (class(phres) != "try-error") {
        if (is.finite(phres$llf)) {
          if (maxllf < phres$llf) {
            maxllf <- phres$llf
            maxph <- ph
           if (verbose$cf1init) cat("o")
          }
          else {
            if (verbose$cf1init) cat("x")
          }
        }
        else {
          if (verbose$cf1init) cat("-")
        }
      }
      else {
        if (verbose$cf1init) cat("-")
      }
    }
    if (verbose$cf1init) cat("\n")
  }
##  if (verbose$cf1init) cat("done\n")
  maxph
}

cf1.param.power <- function(size, mean, s, class) {
  alpha <- rep(1.0/size, size)
  rate <- numeric(size)

  p <- exp(1.0/(size-1) * log(s))
  total <- 1.0
  tmp <- 1.0
  for (i in 2:size) {
    tmp <- tmp * i / ((i-1) * p)
    total <- total + tmp
  }
  base <- total / (size * mean)
  tmp <- base
  for (i in 1:size) {
    rate[i] <- tmp
    tmp <- tmp * p
  }
  cf1(alpha=alpha, rate=rate, class=class)
}

cf1.param.linear <- function(size, mean, s, class) {
  alpha <- rep(1.0/size, size)
  rate <- numeric(size)

  al <- (s-1)/(size-1)
  total <- 1.0
  for (i in 2:size) {
    total <- total + i / (al * (i-1) + 1)
  }
  base <- total / (size * mean)
  for (i in 1:size) {
    tmp <- base * (al * (i-1) + 1)
    rate[i] <- tmp
  }
  cf1(alpha=alpha, rate=rate, class=class)
}

setMethod("emfit.print", signature(model = "cf1"),
  function(model, ...) {
    cat(gettextf("Size : %d\n", model@size))
    cat("Initial : ", model@alpha, "\n")
    cat("Rate    : ", model@rate, "\n")
  }
)

setMethod("emfit.init", signature(model = "cf1"),
  function(model, data, verbose=list(), ...) {
    cf1.param(size=model@size, data=data, verbose=verbose, class=class(model@Q), ...)
  }
)

#### mstep

setMethod("emfit.mstep", signature(model = "cf1"),
  function(model, eres, data, ...) {
    res <- .Call(phfit_mstep_cf1, model, eres, data)
    model@alpha <- res[[1]]
    model@xi <- res[[2]]
    model@Q@x <- res[[3]]
    model@rate <- -diag(model@Q)
    model
  })

