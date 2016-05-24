donostah <- function(x, control)
{
    if(missing(control))
        control <- covRob.control(estim="donostah")
        
    n <- nrow(x)
    p <- ncol(x)

    center <- control$center
    nresamp <- control$nresamp
    maxres <- control$maxres
    prob <- control$prob
    eps <- control$eps

    if(!control$random.sample) 
    {
        if(exists(".Random.seed", where = 1)) 
        {
            random.seed <- get(".Random.seed", pos = 1)
            on.exit(assign(".Random.seed", random.seed, pos = 1))
        }
        set.seed(21)
    }

    if(casefold(nresamp) == "auto")
        nresamp <- ceiling(log(1 - control$prob)/log(1 - (1 - control$eps)^(p+1)))
    else if(!is.integer(nresamp)) 
        stop("nresamp must be a nonnegative integer or ", dQuote("auto"))

    if(nresamp != 0)
        nresamp <- max(1000, nresamp)

    if(casefold(maxres) == "auto")
        maxres <- 2 * nresamp
    else if(!is.integer(maxres))
        stop(sQuote("maxres"), " is not a positive integer")

    tune <- sqrt(qchisq(control$tune, p))

    icent <- 1
    locat <- double(p)
    covmat <- matrix(0.0, p, p)
    storage.mode(covmat) <- "double"
    wk <- double(4*n+p)
    iwork <- integer(4*n+p)
    nresper <- 0
    w <- double(n)
    z <- double(n)

    if(length(center) == 1 && !center)
        center <- rep(0, p)

    if(length(center) > 1) 
    {
        if(length(center) != p)
            stop("the dimension of ", sQuote("center"), " does not match the ", "dimension of ", sQuote("x"))
        x <- sweep(x, 2, center)
        icent <- 0
    }

    sdlist <- .Fortran("rlds",
                      n = as.integer(n),
                      p = as.integer(p),
                      nresamp = as.integer(nresamp),
                      x = as.double(x),
                      tune = as.double(tune),
                      wk = as.double(wk),
                      center = as.double(locat),
                      cov = covmat,
                      maxres = as.integer(maxres),
                      nresper = as.integer(nresper),
                      weights = as.double(w),
                      outlyingness = as.double(z),
                      icent = as.integer(icent),
                      iwork = as.integer(iwork),
                      PACKAGE = "rrcov")

    dist <- mahalanobis(x,
        center = if(length(center) > 1) rep(0, p) else sdlist$center,
        cov = sdlist$cov)

  consistency.correction <- median(dist) / qchisq(.5, p)
    sdlist$cov <- sdlist$cov * consistency.correction
    sdlist$dist <- dist / consistency.correction

    if(length(center) > 1)
        sdlist$center <- center

    sdlist[c("cov", "center", "dist")]
}

covRob.control <- function(estim, ...)
{
    estim <- casefold(estim)
  control <- list(...)
  control$estim <- estim

    if(estim == "donostah") {

    if(is.null(control$nresamp))
      control$nresamp <- "auto"

    if(is.null(control$maxres))
      control$maxres <- "auto"

    if(is.null(control$random.sample))
      control$random.sample <- FALSE

    if(is.null(control$center))
      control$center <- TRUE

    if(is.null(control$tune))
      control$tune <- 0.95

    if(is.null(control$prob))
      control$prob <- 0.99

    if(is.null(control$eps))
      control$eps <- 0.5

    control <- control[c("estim", "nresamp", "maxres", "random.sample",
                         "center", "tune", "prob", "eps")]
  }

    else if(estim == "mcd" || estim == "weighted") {

  ## For backwards compatibility we support the use of quan and ntrial   ##
  ## to specify alpha and nsamp for estim = "mcd", estim = "weighted"    ##
  ## and estim = "M". Providing both quan and alpha or both ntrial and   ##
  ## nsamp will result in an error.                                      ##

    if(is.null(control$alpha))
      control$alpha <- ifelse(is.null(control$quan), 0.5, control$quan)

    if(is.null(control$nsamp))
      control$nsamp <- ifelse(is.null(control$ntrial), 500, control$ntrial)

    if(is.null(control$trace))
      control$trace <- FALSE

    if(is.null(control$use.correction))
      control$use.correction <- TRUE

    if(is.null(control$tolSolve))
      control$tolSolve <- 1e-14

    if(is.null(control$seed))
      control <- control[c("estim", "alpha", "nsamp", "trace",
                           "use.correction", "tolSolve")]
    else
      control <- control[c("estim", "alpha", "nsamp", "seed", "trace",
                           "use.correction", "tolSolve")]
  }

    else if(estim == "m") {

    if(is.null(control$alpha))
      control$alpha <- ifelse(is.null(control$quan), 0.5, control$quan)

    if(is.null(control$nsamp))
      control$nsamp <- ifelse(is.null(control$ntrial), 500, control$ntrial)

    if(is.null(control$trace))
      control$trace <- FALSE

    if(is.null(control$use.correction))
      control$use.correction <- TRUE

    if(is.null(control$tolSolve))
      control$tolSolve <- 1e-14

    if(is.null(control$seed))
      init.control <- control[c("estim", "alpha", "nsamp", "trace",
                                "use.correction", "tolSolve")]
    else
      init.control <- control[c("estim", "alpha", "nsamp", "seed", "trace",
                                "use.correction", "tolSolve")]

    init.control$estim = "mcd"
    control$init.control <- init.control

    if(is.null(control$r))
      control$r <- 0.45

    if(is.null(control$arp))
      control$arp <- 0.05

    if(is.null(control$eps))
      control$eps <- 1e-03

    if(is.null(control$maxiter))
      control$maxiter <- 120

        control <- control[c("estim", "r", "arp", "eps", "maxiter",
                          "init.control")]
  }

    else
        control <- control["estim"]

    control
}
