ramps.control <- function(iter = 1000, beta, sigma2.e, phi, sigma2.z, sigma2.re,
                          z.monitor = TRUE, mpdfun = c("mpdbeta", "mpdbetaz"),
                          file)
{
  iter <- as.integer(iter)
  if (any(iter <= 0)) stop("Only positive integers allowed for 'iter'")
  iter <- if (length(iter) == 1) 1:iter
          else sort(unique(iter))

  if (missing(beta)) beta <- param(NULL)
  else if (!is.param(beta)) stop("Incompatible data type for 'beta'")
  else if (beta$prior != "flat") stop("Only flat priors allowed for 'beta'")

  if (missing(sigma2.e)) sigma2.e <- param(NULL)
  else if (!is.param(sigma2.e)) stop("Incompatible data type for 'sigma2.e'")
  else if (sigma2.e$prior != "invgamma")
    stop("Only inverse gamma priors allowed for 'sigma2.e'")
  if (is.null(sigma2.e$tuning))
    sigma2.e$tuning <- rep(1, length.out = length(sigma2.e))
  else if (any(sigma2.e$tuning <= 0))
    stop("sigma2.e tuning parameters must be > 0")

  if (missing(phi)) phi <- param(NULL)
  else if (!is.param(phi)) stop("Incompatible data type for 'phi'")
  else if (phi$prior == "uniform") {
    if (is.null(phi$tuning)) phi$tuning <- phi$max - phi$min
    else if (any(phi$tuning <= 0)) stop("phi tuning parameters must be > 0")
    else phi$tuning <- (phi$max - phi$min) * phi$tuning
    phi$f <- function(x) prod(dunif(x, phi$min, phi$max))
  }
  else if (phi$prior == "user") {
    if (is.null(phi$tuning))
      stop("phi tuning parameters must supplied")
    else if (any(phi$tuning <= 0))
      stop("phi tuning parameters must be > 0")
  }
  else stop("Only uniform or user-specified priors allowed for 'phi'")

  if (missing(sigma2.z)) sigma2.z <- param(NULL)
  else if (!is.param(sigma2.z)) stop("Incompatible data type for 'sigma2.z'")
  else if (sigma2.z$prior != "invgamma")
    stop("Only inverse gamma priors allowed for 'sigma2.z'")
  if (is.null(sigma2.z$tuning))
    sigma2.z$tuning <- rep(1, length.out = length(sigma2.z))
  else if (any(sigma2.z$tuning <= 0))
    stop("sigma2.z tuning parameters must be > 0")

  if (missing(sigma2.re)) sigma2.re <- param(NULL)
  else if (!is.param(sigma2.re)) stop("Incompatible data type for 'sigma2.re'")
  else if (sigma2.re$prior != "invgamma")
    stop("Only inverse gamma priors allowed for 'sigma2.re'")
  if (is.null(sigma2.re$tuning))
    sigma2.re$tuning <- rep(1, length.out = length(sigma2.re))
  else if (any(sigma2.re$tuning <= 0))
    stop("sigma2.re tuning parameters must be > 0")

  if (!(class(z.monitor) %in% c("logical", "matrix", "data.frame")))
    stop("'z.monitor' must be a logical vector or a dataset of coordinates")

  fnames <- list(params = NULL, z = NULL)
  if (!missing(file)) {
    if (is.list(file)) {
      x <- c("params", "z")
      fnames <- file[match(x, names(file))]
      names(fnames) <- x
    } else {
      file <- file[1:2]
      if (is.null(names(file))) names(file) <- c("params", "z")
      if (!is.na(x <- file["params"])) fnames$params <- as.character(x)
      if (!is.na(x <- file["z"])) fnames$z <- as.character(x)
    }
  }

  list(iter = iter, beta = beta, sigma2.e = sigma2.e, phi = phi,
       sigma2.z = sigma2.z, sigma2.re = sigma2.re,
       z = list(monitor = z.monitor), mpdfun = match.arg(mpdfun),
       file = fnames, expand = 0)
}


param <- function(init, prior = c("flat", "invgamma", "normal", "uniform", "user"),
                  tuning, ...)
{
  retval <- list(init = NULL)
  retval$prior <- match.arg(prior)
  hyper <- list(...)

  n <- length(init)
  switch(retval$prior,
    flat = {
    },
    invgamma = {
      if (!is.numeric(hyper$shape) || hyper$shape <= 0)
        stop("Inverse gamma shape hyperparameter must be numeric > 0")
      if (!is.numeric(hyper$scale) || hyper$scale <= 0)
        stop("Inverse gamma scale hyperparameter must be numeric > 0")
      retval$shape <- rep(hyper$shape, length.out = n)
      retval$scale <- rep(hyper$scale, length.out = n)
      if (any(na <- is.na(init)))
        init[na] <- 1 / rgamma(sum(na), retval$shape[na] / retval$scale[na], 1)
      if (any(init <= 0)) stop("Initial values must be > 0")
    },
    normal = {
      if (!is.numeric(hyper$mean))
        stop("Normal mean hyperparameter must be numeric")
      if (!is.numeric(hyper$var))
        stop("Normal variance hyperparameter must be numeric")
      retval$mean <- rep(hyper$mean, length.out = n)
      if (is.matrix(hyper$var)) {
        val <- hyper$var
      } else {
        val <- diag(n)
        diag(val) <- rep(hyper$var, length.out = n)
      }
      if (any(dim(val) != n)) stop("Non-comformable normal variance matrix")
         val <- try(chol(val))
      if (class(val) == "try-error")
        stop("Normal variance hyperparameter must be positive definite")
      retval$precision <- chol2inv(val)
      if (any(na <- is.na(init)))
        init[na] <- (retval$mean + val %*% rnorm(n))[na]
    },
    uniform = {
      if (!is.numeric(hyper$min))
        stop("Uniform min hyperparameter must be numeric")
      if (!is.numeric(hyper$max))
        stop("Uniform max hyperparameter must be numeric")
      if (any(hyper$min >= hyper$max))
        stop("Uniform min hyperparameter must be < max")
      retval$min <- rep(hyper$min, length.out = n)
      retval$max <- rep(hyper$max, length.out = n)
      if (any(na <- is.na(init)))
        init[na] <- runif(sum(na), retval$min[na], retval$max[na])
      if (any(init <= retval$min, init >= retval$max))
        stop("Initial values must be contained in (min, max)")
    },
    user = {
      if (!is.function(hyper$f))
        stop("User-specified prior function f must be supplied")
      retval$f <- hyper$f
      val <- try(retval$f(init))
      if (length(val) != 1 || !is.finite(val) || val <= 0)
        stop("Prior function f(init) must evaluate to a positive number")
    }
  )

  if (!all(is.finite(init))) stop("Initial values must be numeric")
  retval$init <- init

  if (!missing(tuning)) {
    if (!all(is.finite(tuning))) stop("Tuning values must be numeric")
    retval$tuning <- rep(tuning, length.out = n)
  }

  structure(retval, class = "param")
}


is.param <- function(x)
{
  class(x) == "param"
}


length.param <- function(x)
{
  length(x$init)
}
