
datagen.stsm <- function(n, model = list(), SigmaEV, labels, n0 = 20, freq = 1, 
  old.version = FALSE)
# matrices in 'model' are not checked
{
  # storage matrices

  n <- n + n0
  y <- rep(NA, n)
  a <- matrix(nrow = n, ncol = ncol(model$Z))

  # disturbance terms, innovations

  if (is.null(model$H))
  {
    eps1 <- rep(0, n)
  } else
    eps1 <- rnorm(n, sd = sqrt(model$H[1]))

  if (old.version)
  {
    #requires package 'mvtnorm'
    Meps <- mvtnorm::rmvnorm(n, mean = rep(0, nrow(model$Q)), 
      sigma = model$Q, method = "eigen", pre0.9_9994 = TRUE)
  } else
  {
    #based on MASS::mvrnorm()
    #mvtnorm::rmvnorm uses an expression different from that inferred from theory
    if (missing(SigmaEV))
      SigmaEV <- eigen(model$Q)
    eps <- rnorm(n * ncol(SigmaEV$vectors))
    Meps <- matrix(eps, ncol = n, byrow = TRUE)    
    Meps <- SigmaEV$vectors %*% diag(sqrt(SigmaEV$values)) %*% Meps
    Meps <- t(Meps)
  }

  # iterations to generate data

  if (!is.null(model$a0)) {
    a0 <- model$a0
  } else 
    a0 <- rep(0, ncol(model$Z))
  for (i in seq(n))
  {
    a0 <- a[i,] <- model$T %*% a0 + Meps[i,]
    y[i] <- model$Z %*% a[i,] + eps1[i]
  }

  # output

  if (n0 > 0)
  {
    y <- y[-seq(n0)]
    a <- a[-seq(n0),]
  }

  y <- ts(y, frequency = freq)
  a <- ts(a, frequency = freq)
  if (is.matrix(a))
    a <- a[,diag(model$Q) != 0]

  if (!missing(labels))
  {
    if (length(labels) == ncol(a)) {
      colnames(a) <- labels
    } else
      warning("The length of 'labels' is not equal to the number of", 
        "  components with non-zero variance.")
  }

  list(data = y, components = a)
}
