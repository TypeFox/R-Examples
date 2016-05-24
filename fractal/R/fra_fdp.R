######################################################
## FRACTAL fractionally differenced process functions
##
##   FDSimulate
##   FDWhittle
##
######################################################

###
# FDSimulate
##

"FDSimulate" <- function(delta, innovations.var=1, method="ce", seed=0)
{
  # define local functions
  "FDCirculantEmbedding" <- function(delta, innovations.var=1, n.sample=16, seed=0)
  {
    # check input parameters
    checkScalarType(delta,"numeric")
    checkScalarType(innovations.var,"numeric")
    checkScalarType(n.sample,"integer")
    checkScalarType(seed,"integer")
    seed <- as.integer(abs(seed))

    if (n.sample < 2)
      stop("n.sample must be an integer greater than unity")

    # initialize variables
    ncumsum <- max((2 * delta - 1) %/% 2 + 1, 0)
    delta   <- delta - ncumsum

    # form fdp avcs
    fdp.acvs <- lmACF(lmModel("fdp",delta=delta, innovations.var=innovations.var),
       lag.max=n.sample, type="covariance")@data

    # concatenate reversed fdp acvs
    fdp.acvs <- c(fdp.acvs, rev(fdp.acvs[2:n.sample]))

    # take DFT
    S <- Re(fft(fdp.acvs, inverse=FALSE))

    # use CE method create bootstrap surrogate
    z <- as.vector(itCall("RS_fractal_bootstrap_circulant_embedding",
      S, seed))[seq(n.sample)]
      #COPY=rep(FALSE,2), #CLASSES=c("matrix","integer"),
      #PACKAGE="ifultools"))[seq(n.sample)]

    if (!ncumsum)
      return(z)

    for (i in seq(ncumsum))
      z <- cumsum(z)

    z
  }

  # check input argument types and lengths
  checkVectorType(delta,"numeric")
  checkScalarType(method,"character")
  if (!is.numeric(innovations.var))
    stop("innovations.var must be a numeric object")
  checkScalarType(seed,"integer")
  seed <- as.integer(abs(seed))

  # initialize variables
  n.innov <- length(innovations.var)
  n.delta <- length(delta)

  # perform length check
  if (n.delta > n.innov)
    innovations.var <- c(innovations.var, rep(innovations.var[n.innov],
      n.delta - n.innov))

  # check method
  method <- match.arg(method, c("ce","cholesky"))

  #TODO: axe this probably: we may only support ce method
  if (method == "cholesky"){
    z <- as.vector(itCall("RS_fractal_fdp_simulate", delta, innovations.var))
      #COPY = rep(FALSE,2), #CLASSES= c("numeric", "numeric"),
      #PACKAGE="ifultools"))
  } else {

    # find unique delta levels and match against
    # original delta vector
    delta.levels <- sort(unique(delta))
    idelta <- match(delta, delta.levels)

    # create matrix of FD realizations, each column
    # corresponds to a unique delta
    z <- apply(matrix(delta.levels), MARGIN=1, function(x, n.delta, seed, FDCirculantEmbedding)
      FDCirculantEmbedding(delta=x, innovations.var=1, n.sample=n.delta),
      n.delta=n.delta, seed=seed, FDCirculantEmbedding=FDCirculantEmbedding)

    z <- z[cbind(seq(n.delta), idelta)] * innovations.var
  }

  oldClass(z) <- "FDSimulate"
  attr(z, "delta")  <- delta
  attr(z, "innov")  <- innovations.var
  attr(z, "method") <- "circulant embedding"

  z
}

##
# FDWhittle
##

"FDWhittle" <- function(x, method="continuous", dc=FALSE,
  freq.max=0.5, delta.min=-1,delta.max=2.5, sdf.method="direct", ...)
{
  # define local functions
  Isumd <- function(delta, Sx, discrete)
  {
  	# function to evaluate the weighted sum of an FD process spectrum
  	f  <- Sx$frequency
  	N  <- Sx$n.sample
  	Np <- (N - 1) %/% 2
    Q  <- 2 * mean(Sx$data * (4 * (sin(pi * f))^2)^delta)
    ifelse1(!discrete, Q, log(Q) + 2 * delta * (((N/(2 * Np)) - 1) *
      log(2) - log(2 * N)/(2 * Np)))
  }

  # check input arguments
  checkScalarType(method,"character")
  method <- match.arg(lowerCase(method),c("discrete","continuous"))
  checkScalarType(freq.max,"numeric")
  if (freq.max <= 0.0 || freq.max > 0.5)
    stop("freq.max must be on the normalized frequency range (0,0.5)")

  # calculate single-sided SDF estimate with a normalized
  # spectral resolution of 1/N
  sdf <- sapa::SDF(x, method=sdf.method, single.sided=TRUE, npad=length(x), ...)

  # associate freq.max normalized frequency
  # with corresponding index into SDF
  # frequency vector
  xatt          <- attributes(sdf)
  freq.indx.max <- ceiling(freq.max/xatt$deltaf)
  freq.indx.min <- ifelse1(dc, 1, 2)
  ifreq         <- seq(freq.indx.min, freq.indx.max)

  # pack relevant spectral information into a list.
  # limit the range of sdf and corresponding frequency
  # according to the input specifications
  Sx <- list(data=sdf[ifreq], frequency=xatt$frequency[ifreq], n.sample=xatt$n.sample)

  # use optimize() function to find the best value of delta
  # return Hurst coefficient estimate
  #dfinal$minimum + 0.5
  z <- optimize(Isumd, lower=delta.min, upper=delta.max, tol=0.001,
    maximum=FALSE, Sx=Sx, discrete=is.element(method,"discrete"))$minimum

  names(z) <- "delta"
  z
}
