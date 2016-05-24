# Compute the usual unbiased estimate of the variance in a linear model.
estimate_noise <- function(X, y, intercept = TRUE) {
  n = nrow(X)
  if (intercept)
    X = cbind(rep(1,n), X)
  p = ncol(X)
  stopifnot(n > p)

  fit = lm.fit(X, y)
  sqrt(sum(fit$residuals^2) / (n-p))
}

# A variant of the built-in function 'scale' that using L2-normalization for
# the scaling.
normalize <- function(X, center = TRUE, scale = TRUE) {
  X = as.matrix(X)
  if (center) {
    means = colMeans(X, na.rm = TRUE)
    X = sweep(X, 2, means)
  }
  if (scale) {
    scales = apply(X, 2, function (x) sqrt(sum(x^2)))
    X = sweep(X, 2, scales, "/")
  }
  X
}

# Generate a random synthetic model and data. Used for testing.
random_problem <- function(n, p, k=NULL, amplitude=3, sigma=1) {
  if (is.null(k))
    k = max(1, as.integer(p/5))

  X = matrix(rnorm(n*p), n, p)
  nonzero = sample(p, k)
  beta = amplitude * (1:p %in% nonzero)
  y = X %*% beta + rnorm(n, sd = sigma)
  list(X = X, y = y, beta = beta)
}

# Evaluate an expression with the given random seed, then restore the old seed.
with_seed <- function(seed, expr) {
  seed.old <- if (exists('.Random.seed')) .Random.seed else NULL
  set.seed(seed)
  on.exit({
    if (is.null(seed.old)) {
      if (exists('.Random.seed'))
        rm(.Random.seed, envir=.GlobalEnv)
    } else {
      .Random.seed <<- seed.old
    }
  })
  expr
}