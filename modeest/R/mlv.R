# Author: P. PONCET

#! a voir : package 'ks', notamment pour les kernel density derivative estimates

.methodList <- c("asselin",
                 #"chernoff",
                 "density",
                 "discrete",
                 #"devroye",
                 #"ekblom",
                 "eme",
                 "grenander",
                 "hrm",
                 "hsm",
                 "kernel",
                 #"kim",
                 "lientz",
                 #"lms",
                 #"logspline", # cf. package logspline
                 "meanshift",
                 "mfv",# most frequent value
                 #"mizoguchi",
                 "naive",
                 "parzen",
                 #"robertson",
                 "shorth",
                 #"splines",                 
                 "tsybakov",
                 "venter", 
                 "vieu")
                 #"wavelet",
                 #"wtp")


# Computes an estimate of the mode of a univariate (unimodal ?) distribution
mlv <-
function(x,
         ...)
{
  UseMethod("mlv")
}


# S3 method for class 'character'
mlv.character <-
function(x,
         ...)
{
  if (!is.character(x)) stop("argument 'x' must be a character")
  x <- tolower(x)
  if (x %in% c("chisquare", "chisquared")) x <- "chisq"
  if (x == "exponential") x <- "exp"
  if (x %in% c("generalised_hyperbolic", "generalized_hyperbolic")) x <- "gh"
  if (x == "fdist") x <- "f"
  if (x == "gammadist") x <- "gamma"
  if (x %in% c("gaussian", "normal")) x <- "norm"
  if (x %in% c("generalised_extreme_value", "generalized_extreme_value")) x <- "gev"
  if (x %in% c("generalised_pareto", "generalized_pareto")) x <- "gpd"
  if (x %in% c("lognormal", "loggaussian")) x <- "lnorm"
  if (x == "hyperbolic") x <- "hyp"
  if (x == "kumaraswamy") x <- "kumar"
  if (x == "logistic") x <- "logis"
  if (x == "normal_inverse") x <- "nig"
  #if (x == "symmetric_stable") x <- "symstb"
  if (x == "student") x <- "t"
  if (x == "uniform") x <- "unif"
  if (x == "bernoulli") x <- "bern"
  if (x == "binomial") x <- "binom"
  if (x == "geometric") x <- "geom"
  if (x %in% c("hypergeometric", "hypergeom")) x <- "hyper"
  if (x == "negative_binomial") x <- "nbinom"
  if (x == "poisson") x <- "pois"
  x <- match.arg(x, .distribList)  # '.distribList' is defined in 'distribMode.R'
  M <- do.call(paste(".mlv.", x, sep = ""), list(...))
  if (inherits(M, "try-error")) {
    stop(paste("argument 'x' is incorrect. Function .mlv.", x, " does not exist", sep = ""))
  }

  th <- M$M

  ## Output
  return(structure(list(M = th,
                        skewness = M$skewness,
                        x = x,
                        method = M$method,
                        bw = NULL, 
                        boot = FALSE,
                        #boot.M = th,
                        call = match.call()),
                   class = "mlv"))
}


# S3 method for class 'factor'
mlv.factor <-
function(x,               # sample (the data)
         ...)
{

  ## Preliminaries
  if (!is.factor(x)) stop("argument 'x' must be a factor")

  th <- discrete(x, ...)

  ## Bickel's measure of skewness
  xx <- as.numeric(x)
  cdf.value <- (length(xx[xx < th]) + length(xx[xx == th])/2)/length(xx)
  skewness <- 1-2*cdf.value

  ## Output
  return(structure(list(M = th,
                        skewness = skewness,
                        x = x,
                        method = "discrete", 
                        bw = NULL, 
                        boot = FALSE,
                        call = match.call()),
                   class = "mlv"))


}


# S3 method for class 'integer'
mlv.integer <-
function(x,               # sample (the data)
         na.rm = FALSE,   # Should NA values in 'x' be removed?
         ...)
{

  ## Preliminaries
  if (!is.integer(x)) stop("argument 'x' must be an integer")
  x <- as.vector(x)

  x.na <- is.na(x)
  if (any(x.na)) {
    if (na.rm) {
      x <- x[!x.na]
    } else {
      stop("argument 'x' contains missing values")
    }
  }

  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
  }

  th <- discrete(x, ...)

  ## Bickel's measure of skewness
  cdf.value <- (length(x[x < th]) + length(x[x == th])/2)/length(x)
  skewness <- 1-2*cdf.value

  ## Output
  return(structure(list(M = th,
                        skewness = skewness,
                        x = x,
                        method = "discrete",
                        bw = NULL, 
                        boot = FALSE,
                        call = match.call()),
                   class = "mlv"))


}


mlv.default <-
function(x,               # sample (the data)
         bw = NULL,       # bandwidth
         method,          # must belong to '.methodList'
         na.rm = FALSE,   # Should NA values in 'x' be removed?
         #dip.test = FALSE,# Should Hartigan's DIP test of unimodality be done?
         boot = FALSE,    # Should the mode be bootstrapped?
         R = 100,         # How many bootstrap resamplings?
         B = length(x),   # Size of the bootstrap samples drawn from 'x'
         ...)
{

  ## Preliminaries
  if (!is.numeric(x)) stop("argument 'x' must be numeric")
  x <- as.vector(x)

  x.na <- is.na(x)
  if (any(x.na)) {
    if (na.rm) {
      x <- x[!x.na]
    } else {
      stop("argument 'x' contains missing values")
    }
  }

  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
  }

  #if (!dip.test) {
  #  dip.test <- NULL
  #} else {
  #  dip.test <- dip(x, full.result = FALSE)
  #  cat("function 'dip' from package 'diptest' is used to compute Hartigan's DIP statistic\n")
  #}

  if (missing(method)) {
    warning("argument 'method' is missing. Data are supposed to be continuous. Default method 'shorth' is used")
    method <- "shorth"
  } else if (pmatch(tolower(method), c("density", "kernel"), nomatch = 0)) {
    method <- "parzen"
  } else method <- match.arg(tolower(method), .methodList) # '.methodList' is defined in 'mlv.R'
  
  if (method == "lientz") method <- "mlv.lientz"
  
  ## Without bootstrapping
  if (!boot) {
    theta <- eval(parse(text = paste(method, "(x, bw = bw, ...)", sep="")))
                                                                            
  ## With bootstrapping
  } else {
    f <- paste(method, "(b, bw = bw, ...)", sep = "")
    theta <- sapply(1:R, function(z)                        #! VOIR AUSSI le package 'boot' !
                         {
                           b <- sample(x, B, replace = TRUE)
                           eval(parse(text = f))
                         })
  }

  th <- mean(theta)

  ## Bickel's measure of skewness
  cdf.value <- (length(x[x < th]) + length(x[x == th])/2)/length(x)
  skewness <- 1-2*cdf.value

  ## Output
  return(structure(list(M = th,
                        skewness = skewness,
                        x = x,
                        method = method,
                        bw = bw, 
                        #dip.stat = dip.test,
                        boot = boot,
                        boot.M = theta,
                        call = match.call()),
                   class = "mlv"))


}


# S3 method for class 'density' (this function is derived from package hdrcde)
mlv.density <-
function(x,
         all = TRUE, 
         #dip.test = FALSE,
         abc = FALSE,
         ...)
{
  if (!inherits(x, "density")) stop("argument 'x' must inherit from class 'density'")
  
  y <- x$y
  x <- x$x

  #if (!dip.test) {
  #  dip.test <- NULL
  #} else {
  #  #require(diptest)
  #  dip.test <- dip(x, full.result = FALSE)
  #}
  
  idx <- y == max(y)
  M <- x[idx]

  ## Bickel's measure of skewness
  cdf.value <- (length(x[x < M]) + length(x[x == M])/2)/length(x)
  skewness <- 1-2*cdf.value
  
  if (all) {
    yy <- c(0, y, 0)
    ny <- length(yy)
    idx <- (yy[2:(ny - 1)] > yy[1:(ny - 2)]) & (yy[2:(ny - 1)] > yy[3:ny])
    M <- unique(x[idx], M)
  }
   
  ## Output
  return(structure(list(M = M,
                        skewness = skewness,
                        x = x,
                        method = "density",
                        #dip.stat = dip.test,
                        boot = FALSE,
                        call = match.call()),
                   class = "mlv"))
}


plot.mlv <-
function(x, # an object of class 'mlv'
         ...)
{
  if (!inherits(x, "mlv")) stop("argument 'x' must inherit from class 'mlv'")
  th <- x$M
  xx <- x$x
  method <- x$method

  arg <- list(...)
  main <- arg$main
  xlab <- arg$xlab
  ylab <- arg$ylab

  if (is.null(main)) main <- deparse(x$call)
  if (is.null(xlab)) xlab <- paste("Mode = ", th, sep="")
  if (is.null(ylab)) ylab <- "Density and mode"

  if (x$boot) {
    layout(matrix(1:2, 1, 2))  
  }
  
  if (method == "discrete") {
    if (is.character(xx)) {
      #p <- paste("curve(d", xx, "(x))") 
      #! à faire
    } else {
      hist.default(xx, freq = TRUE, main = main, xlab = xlab, ylab = ylab, ...)    
      abline(v = th, col = 2)
    }

  } else {
    if (is.character(xx)) {
      #p <- paste("curve(d", xx, "(x))")
      #! à faire
    } else {
      plot.default(density.default(xx), main = main, xlab = xlab, ylab = ylab, ...)
      abline(v = th, col = 2)
    }    
  }
  
  if (x$boot) {
    boxplot.default(x$boot.M, main = "Bootstrap of the mode: related boxplot")
  }
  
  layout(matrix(1, 1, 1))  
  return(invisible(NULL))
}


print.mlv <-
function(x, # an object of class 'mlv'
         digits = NULL,
         ...)
{
  if (!inherits(x, "mlv")) stop("argument 'x' must inherit from class 'mlv'")
  th <- x$M
  #xx <- x$x
  skew <- x$skewness
  method <- x$method
  #dip.stat <- x$dip.stat

  if (method == "discrete") {
    #cat("Nature of the data: discrete\n")
    cat("Mode (most frequent value):", format(th, digits = digits), "\n")
    cat("Bickel's modal skewness:", format(skew, digits = digits), "\n")

  } else {
    #cat("Nature of the data: continuous\n")
    cat("Mode (most likely value):", format(th, digits = digits), "\n")
    cat("Bickel's modal skewness:", format(skew, digits = digits), "\n")
  }
  
  #if (!is.null(dip.stat)) cat("Hartigan's DIP statistic:", format(dip.stat, digits = digits), "\n")
  cat("Call:", deparse(x$call), "\n")

  return(invisible(th))
}

