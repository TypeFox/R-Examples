#***********************************************************************
## Transformation of a vector of "Renouv" parameters to Gumebl block max
## and reciprocal transform.
##
## Author: Yves Deville <deville.yves@alpestat.com> 
##
##***********************************************************************

##=======================================================================
## Here 'w' can be a vector, in which case the Gumbel parameters
## are returned as a matrix.
##=======================================================================

Ren2gumbel <- function(object,
                       threshold = NULL,
                       w = 1,
                       distname.y = c("exponential", "exp"),
                       jacobian = (length(w) == 1L),
                       vcovRen = NULL) {
  
  if (inherits(object, "Renouv")) {
    threshold <- object$threshold
    distname.y <- object$distname.y
    if (!distname.y %in% c("exp", "exponential")) {
      stop("'distname.y' must be \"exp\" or \"exponential\"")
    }
    parRen <- coef(object)
    vcovRen <- vcov(object)
    jacobian <- TRUE
  } else {
    parRen <- object
    if (is.null(threshold)) stop("'threshold' must be given")
    distname.y <- match.arg(distname.y)
  }
  
  RenNames <- names(parRen)
  
  ## now convert to Gumbel
  parGumbel <- array(NA, dim = c(length(w), 2L),
                  dimnames = list(NULL, c("loc", "scale")))

  lambda <- parRen["lambda"]
  LW <- lambda * w
  logLW <- log(LW)
  parGumbel[ , "loc"] <- threshold  +  logLW / parRen["rate"] 
  parGumbel[ , "scale"] <- 1.0 / parRen["rate"]  
  
  if (jacobian) {
    
    if (length(w) > 1L)
      stop("jacobian computation is only possible when length(w) is 1")
   
    ## The jacobian
    derParGumbel <- matrix(0, nrow = 2L, ncol = 2L)
    rownames(derParGumbel) <- c("loc", "scale")
    colnames(derParGumbel) <- c("lambda", "rate")
    ## row 1
    derParGumbel["loc", "lambda"] <- 1.0 / lambda /  parRen["rate"] 
    derParGumbel["loc", "rate"] <-  -logLW / parRen["rate"] / parRen["rate"]
    ## row 2
    derParGumbel["scale", "lambda"] <- 0.0 
    derParGumbel["scale", "rate"] <- - 1.0 / parRen["rate"] / parRen["rate"]
    
    attr(parGumbel, "jacobian") <- derParGumbel
    
  }
  
  ## covariance matrix
  if (!is.null(vcovRen)) {

    if (!jacobian) stop("when 'vcovRen' is provided 'jacobian' must be TRUE")
    
    if ( !is.matrix(vcovRen) || (nrow(vcovRen) != 2L) || (ncol(vcovRen) != 2L) ) {
      stop("'vcovRen' must be a matrix with 2 rows and 2 columns")
    }
    if ( !setequal(rownames(vcovRen), RenNames) ) {
      stop("'vcovRen' must have suitable rownames and colnames, see help")
    }
    
    colnames(vcovRen) <- rownames(vcovRen)
    vcovGumbel <- derParGumbel %*% vcovRen[RenNames, RenNames] %*% t(derParGumbel)
    attr(parGumbel, "vcov") <- vcovGumbel
  }
  
  if (length(w) == 1L) {
    parGumbel <- drop(parGumbel)
  }
  
  attr(parGumbel, "threshold") <- threshold
  attr(parGumbel, "jacobian") <- derParGumbel

  parGumbel
  
}

##=======================================================================
## Here 'w' must be of length 1
##=======================================================================

gumbel2Ren <- function(parGumbel,
                       threshold = NULL,
                       lambda = NULL,
                       w = 1,
                       distname.y = c("exponential"),
                       vcovGumbel = NULL,
                       jacobian = TRUE,
                       plot = FALSE) {
  
  gumbelNames <- c("loc", "scale")
  if ( !setequal(names(parGumbel), gumbelNames) ) {
    stop("'parGumbel' must have suitable names, see help")
  }
  
  distname.y <- match.arg(distname.y)
  parRen <- c("lambda" = NA, "rate" = NA)

  if (!is.null(threshold)) {

    fixed <- "threshold"

    if (!is.null(lambda)) {
      stop("only one of 'threshold' and 'lambda' can be given")
    }
    
    C <- ( threshold - parGumbel["loc"] ) / parGumbel["scale"]
    ## if (C <= 0) stop("'threshold' too small. Must be >= ", parGumbel["loc"])
   
    E <- exp(-C)
    lambda <- E / w
    LW <- lambda * w
    logLW <- log(LW)

  } else {
    fixed <- "lambda"
    ## 'lambda' provided
    LW <- lambda * w
    logLW <- log(LW)
    threshold <-  parGumbel["loc"] - logLW * parGumbel["scale"]
  }
  
  parRen["lambda"] <- lambda
  parRen["rate"] <- 1.0 / parGumbel["scale"]

  ## the plot uses untransformed parameters.

  if (plot) {
    
    lims <- qgumbel(c(0.01, 0.9999), loc = parGumbel["loc"],
                 scale = parGumbel["scale"])
    
    xlev <- seq(from = lims[1], to = lims[2], length.out = 1000)

    T.gumbel <- 1 / pgumbel(xlev, loc = parGumbel["loc"], scale = parGumbel["scale"],
                            lower.tail = FALSE)

    ## Caution: no "loc" parameter here! 
    T.exp <- 1 / parRen["lambda"] /
      pexp(xlev - threshold, rate = parRen["rate"], lower.tail = FALSE)
    
    ## exponential return level plot
    plot(log(T.gumbel, 10), xlev, type = "l",  lty = "solid")
    lines(log(T.exp, 10), xlev, type = "l", lty = "dashed", col = "red")
    
  }
  
  if (jacobian) {
    
    derParRen  <- matrix(0, nrow = 3,  ncol = 2)
    colnames(derParRen) <- c("loc", "scale")
    rownames(derParRen) <- c("lambda", "threshold", "rate")
    
    if (fixed == "threshold") {      
      ## auxiliary variables 'C' and 'E' defined above 
      derParRen["lambda", ] <- c(1.0, C) *  E / w / parGumbel["scale"]
      derParRen["threshold", "loc"] <- 1.0
      ## derParRen["rate", "scale" ] <- 0.0   ## done 
      derParRen["rate", "scale" ] <- - 1.0 / parGumbel["scale"] / parGumbel["scale"]
      
    } else if (fixed == "lambda") {
      derParRen["threshold", "loc"] <- 1.0
      derParRen["threshold", "scale"] <- - logLW
      derParRen["rate", "loc"] <- 0.0
      derParRen["rate", "scale"] <- - 1.0 / parGumbel["scale"] / parGumbel["scale"]
    }

    attr(parRen, "jacobian") <- derParRen

  }

  if (!is.null(vcovGumbel)) {
    if ( !is.matrix(vcovGumbel) || (nrow(vcovGumbel) != 2L) || (ncol(vcovGumbel) != 2L) ) {
      stop("'vcovGumbel' must be a matrix with 2 rows and 2 columns")
    }
    if ( !setequal(rownames(vcovGumbel), gumbelNames) ) {
      stop("'vcovGumbel' must have suitable rownames and colnames, see help")
    }
    
    colnames(vcovGumbel) <- rownames(vcovGumbel)
    
    print(vcovGumbel[gumbelNames, gumbelNames])
    print(derParRen)
    
    vcovRen <- derParRen %*% vcovGumbel[gumbelNames, gumbelNames] %*% t(derParRen)
    print(vcovRen[c("lambda", "rate"), c("lambda", "rate")])
    
    attr(parRen, "vcov") <- vcovRen
  } 
  
  attr(parRen, "threshold") <- as.numeric(threshold)
  attr(parRen, "distname.y") <- as.character(distname.y)
  
  parRen    
  
}
