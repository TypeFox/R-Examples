rmnorm <- 
function(n = 1, mean = rep(0, nrow(Sigma)), Sigma = diag(length(mean)))
{
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(mean) != nrow(Sigma))
    stop("mean and sigma have non-conforming size")
  p <- nrow(Sigma)

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  # call C code
  y <- .C("rand_norm", 
          y = as.double(y),
          dims = as.integer(dy),
          mean = as.double(mean),
          Sigma = as.double(Sigma))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}

rmCauchy <- 
function(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)))
{
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(center) != nrow(Scatter))
    stop("center and scatter have non-conforming size")
  p <- nrow(Scatter)

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  # call C code
  y <- .C("rand_cauchy", 
          y = as.double(y),
          dims = as.integer(dy),
          center = as.double(center),
          Scatter = as.double(Scatter))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}

rmt <- 
function(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)), df = 4)
{
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(center) != nrow(Scatter))
    stop("center and scatter have non-conforming size")
  if (df <= 0)
    stop("Degrees of freedom must be > 0")
  p <- nrow(Scatter)

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  # call C code
  y <- .C("rand_student", 
          y = as.double(y),
          dims = as.integer(dy),
          center = as.double(center),
          Scatter = as.double(Scatter),
          df = as.double(df))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}

rmslash <- 
function(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)), df = 2)
{
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(center) != nrow(Scatter))
    stop("center and scatter have non-conforming size")
  if (df <= 0)
    stop("Degrees of freedom must be > 0")
  p <- nrow(Scatter)

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  # call C code
  y <- .C("rand_slash", 
          y = as.double(y),
          dims = as.integer(dy),
          center = as.double(center),
          Scatter = as.double(Scatter),
          df = as.double(df))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}

rmcnorm <- 
function(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)), 
  epsilon = 0.05, vif = 0.25)
{
  if (n <= 0)
    stop("n must be a positive integer")
  if (length(center) != nrow(Scatter))
    stop("center and scatter have non-conforming size")
  if ((epsilon < 0) || (epsilon > 1))
      stop("contamination percentage must be in [0,1]")
  if ((vif <= 0) || (vif >= 1))
      stop("variance inflation factor must be in (0,1)")
  p <- nrow(Scatter)

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  # call C code
  y <- .C("rand_contaminated", 
          y = as.double(y),
          dims = as.integer(dy),
          center = as.double(center),
          Scatter = as.double(Scatter),
          epsilon = as.double(epsilon),
          vif = as.double(vif))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}

rsphere <- 
function(n = 1, p = 2)
{
  if (n <= 0)
    stop("n must be a positive integer")
  if (p <= 1)
    stop("dimension must be > 1")

  y <- matrix(0, nrow = n, ncol = p)
  dy <- dim(y)
  # call C code
  y <- .C("rand_sphere", 
          y = as.double(y),
          dims = as.integer(dy))$y
  y <- matrix(y, nrow = n, byrow = TRUE)
  y
}
