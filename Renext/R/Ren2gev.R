##***********************************************************************
## Transformation of a vector of "Renouv" parameters to gev block max
## and reciprocal transform.
##
## Author: Yves Deville <deville.yves@alpestat.com> 
##
##***********************************************************************


##=======================================================================
## Here 'w' can be a vector, in which case the GEV parameters
## are returned as a matrix.
##=======================================================================

Ren2gev <-
function(object,
         threshold = NULL, w = 1,
         distname.y = c("gpd", "GPD", "lomax", "maxlo"),
         jacobian = (length(w) == 1L),
         vcovRen = NULL) {

  if (inherits(object, "Renouv")) {
    threshold <- object$threshold
    distname.y <- object$distname.y
    if (!distname.y %in% c("gpd", "GPD", "lomax", "maxlo")) {
      stop("'distname.y' must be \"gpd\", \"GPD\", \"lomax\" or \"maxlo\"")
    }
    parRen <- coef(object)
    vcovRen <- vcov(object)
    jacobian <- TRUE
  } else {
    parRen <- object
    if (is.null(threshold)) stop("'threshold' must be given")
    distname.y <- match.arg(distname.y)
  }
  
  ## convert to gpd if needed
  if ( !(distname.y %in% c("gpd","GPD")) ) {
    
    if (distname.y == "lomax") {
      transParRen <- lomax2gpd(parRen[-1L], jacobian = TRUE)
    } else if (distname.y == "maxlo") {
      transParRen <- maxlo2gpd(parRen[-1L], jacobian = TRUE)
    } else {
      stop("Only \"gpd\", \"lomax\" and \"maxlo\" are accepted distributions")
    }
    transNames <- names(transParRen)  ## gpd names!
    parRen <- c("lambda" = as.numeric(parRen["lambda"]), transParRen)
    
    ## cat("parRen apres transform\n"); print(parRen)      
    
  }
  
  RenNames <- names(parRen)
  ## cat("RenNames\n")
  ## print(RenNames)
  ## print(vcovRen)
  
  ## now convert to GEV
  parGev <- array(NA, dim = c(length(w), 3L),
                  dimnames = list(NULL, c("loc", "scale", "shape")))
  
  prov <- (parRen["lambda"] * w) ^ parRen["shape"]
  parGev[ , "loc"] <- threshold  + (prov - 1) * parRen["scale"] / parRen["shape"] 
  parGev[ , "scale"] <- parRen["scale"] * prov 
  parGev[ , "shape"] <- parRen["shape"]
  
  if (jacobian) {
    
    if (length(w) > 1L)
      stop("jacobian computation is only possible when length(w) is 1")
    
    L <- parRen["lambda"] * w
    Le <- (parRen["lambda"] * w)^parRen["shape"]
    LBoxCox <- (Le - 1) / parRen["shape"]
    
    ## The jacobian id for GPD parameters
    derParGev <- matrix(0, nrow = 3L, ncol = 3L)
    rownames(derParGev) <- c("loc", "scale", "shape")
    colnames(derParGev) <- c("lambda", "scale", "shape")
    ## row 1
    derParGev["loc", "lambda"] <- Le * parRen["scale"] / parRen["lambda"]
    derParGev["loc", "scale"] <- LBoxCox
    derParGev["loc", "shape"] <-
      (Le *log(L) - LBoxCox) * parRen["scale"] / parRen["shape"] 
    ## row 2
    derParGev["scale", "lambda"] <- Le * parRen["shape"] * parRen["scale"] / parRen["lambda"]
    derParGev["scale", "scale"] <- Le
    derParGev["scale", "shape"] <- Le *log(L) * parRen["scale"]
    ## row 3
    ## derParGev["shape", "lambda"] <- 0.0       ### remind...
    ## derParGev["shape", "scale"] <- 0.0        ### remind... 
    derParGev["shape", "shape"] <- 1.0

    if ( !(distname.y %in% c("gpd", "GPD")) ) {
      derTrans  <- matrix(0, nrow = 3L,  ncol = 3L)
      colnames(derTrans) <- c("lambda", transNames)
      rownames(derTrans) <- c("lambda", "scale", "shape")
      derTrans["lambda", "lambda"] <- 1.0
      derTrans[c("scale", "shape"), transNames] <-
        attr(transParRen, "jacobian")
      derParGev <- derParGev %*%  derTrans 
    }
    
    attr(parGev, "jacobian") <- derParGev
    
  }
  
  ## covariance matrix

  if (!is.null(vcovRen)) {

    if (!jacobian) stop("when 'vcovRen' is provided 'jacobian' must be TRUE")
    
    if ( !is.matrix(vcovRen) || (nrow(vcovRen) != 3L) || (ncol(vcovRen) != 3L) ) {
      stop("'vcovRen' must be a matrix with 3 rows and 3 columns")
    }
    if ( !setequal(rownames(vcovRen), RenNames) ) {
      stop("'vcovRen' must have suitable rownames and colnames, see help")
    }
    
    colnames(vcovRen) <- rownames(vcovRen)
    vcovGev <- derParGev %*% vcovRen[RenNames, RenNames] %*% t(derParGev)
    attr(parGev, "vcov") <- vcovGev
  }
  
  if (length(w) == 1L) {
    parGev <- drop(parGev)
  }
  
  attr(parGev, "threshold") <- threshold
  attr(parGev, "jacobian") <- derParGev

  parGev
  
}

##=======================================================================
## Here 'w' must be of length 1
##=======================================================================

gev2Ren <-
function(parGev,
         threshold = NULL,
         lambda = NULL,
         w = 1,
         distname.y = c("gpd", "GPD", "lomax", "maxlo"),
         vcovGev = NULL,
         jacobian = TRUE,
         plot = FALSE) {
  
  gevNames <- c("loc", "scale", "shape")
  if ( !setequal(names(parGev), gevNames) ) {
    stop("'parGev' must have suitable names, see help")
  }
  
  distname.y <- match.arg(distname.y)
  if (parGev["shape"] <= 0 && distname.y == "lomax") 
    stop("GEV shape parameter <= 0 can not be used with \"lomax\"")
  
  if (parGev["shape"] >= 0 && distname.y == "maxlo") 
    stop("GEV shape parameter >= 0 can not be used with \"maxlo\"")
  
  parRen <- c("lambda" = NA, "scale" = NA, "shape" = NA)

  if (!is.null(threshold)) {

    fixed <- "threshold"
    ## 'threshold' provided
    umod <- (threshold - parGev["loc"]) / parGev["scale"]
    C <- (1 + parGev["shape"] * umod)
    
    if (C <= 0) {
      lim <- parGev["loc"] - parGev["scale"] / parGev["shape"]
      if (parGev["shape"] > 0) stop("'threshold' too small. Must be >= ", lim)
      else if (parGev["shape"] < 0) stop("'threshold' too large. Must be >= ", lim)
    }
    
    if (!is.null(lambda)) {
      stop("only one of 'threshold' and 'lambda' can be given")
    }
    lambda <- C^(-1/parGev["shape"])  / w
    L <- lambda * w

  } else {
    fixed <- "lambda"
    ## 'lambda' provided
    L <- lambda * w
    Le <- L ^(-parGev["shape"])
    threshold <-  parGev["loc"] +
      parGev["scale"] *( L^(-parGev["shape"]) -1 ) / parGev["shape"]
  }
  
  parRen["lambda"] <- lambda
  parRen["shape"] <- parGev["shape"]
  parRen["scale"] <- parGev["scale"] / L^parGev["shape"]

  ## the plot uses untransformed parameters.

  if (plot) {
    
    lims <- qgev(c(0.01, 0.9999), loc = parGev["loc"],
                 scale = parGev["scale"], shape = parGev["shape"])
    
    xlev <- seq(from = lims[1], to = lims[2], length.out = 1000)

    T.gev <- 1 / pgev(xlev, loc = parGev["loc"], scale = parGev["scale"],
                      shape = parGev["shape"], lower.tail = FALSE)
    
    T.gpd <- 1 / parRen["lambda"] /
      pgpd(xlev, loc = threshold, scale = parRen["scale"], shape = parRen["shape"],
           lower.tail = FALSE)
    
    ## exponential return level plot
    plot(log(T.gev, 10), xlev, type = "l",  lty = "solid")
    lines(log(T.gpd, 10), xlev, type = "l", lty = "dashed", col = "red")
    
  }
  
  ## transform 
  if ( !(distname.y %in% c("gpd", "GPD")) ) {

    if (distname.y == "lomax") {
      transParRen <- gpd2lomax(parRen[c("scale", "shape")], jacobian = TRUE)
    } else if (distname.y == "maxlo") {
      transParRen <- gpd2maxlo(parRen[c("scale", "shape")], jacobian = TRUE)
    }
    transNames <- names(transParRen)
    parRen <- c("lambda" = as.numeric(parRen["lambda"]), transParRen)
    ## cat("XXXX", transNames, "\n")
  }
  
  if (jacobian) {
    
    derParRen  <- matrix(0, nrow = 4,  ncol = 3)
    colnames(derParRen) <- c("loc", "scale", "shape")
    rownames(derParRen) <- c("lambda", "threshold", "scale", "shape")
    ##cat("fixed = ", fixed, "\n")
    
    if (fixed == "threshold") {      
      ## auxiliary variables 'umod' and 'C' defined above 
      derC <- c("loc" = - parGev["shape"] / parGev["scale"],
                "scale" = - umod * parGev["shape"] / parGev["scale"],
                "shape" = umod)
      
      derParRen["lambda", ] <- -C^(-1/parGev["shape"] - 1) * derC / w / parGev["shape"]
      derParRen["lambda", "shape"] <- derParRen["lambda", "shape"] +
        lambda * log(C) / parGev["shape"] / parGev["shape"]
      derParRen["scale", ] <- parGev["scale" ] * derC
      derParRen["scale", "scale" ] <-  derParRen["scale", "scale" ] + C
      derParRen["shape", ] <- c(0, 0, 1)
      
    } else if (fixed == "lambda") {
      ## auxiliary variables 'L' and 'Le' defined above
      ## cat(sprintf("L = %5.2f Le = %5.2f\n", L, Le)) 
      derParRen["threshold", "loc"] <- 1.0
      derParRen["threshold", "scale"] <- -(1 - Le) / parGev["shape"]
      derParRen["threshold", "shape"] <-
        ( (1 - Le) / parGev["shape"] - log(L) * Le ) * parGev["scale"] / parGev["shape"]
      derParRen["scale", "loc"] <- 0.0
      derParRen["scale", "scale"] <- Le
      derParRen["scale", "shape"] <- - log(L) * Le * parGev["scale"]
      derParRen["shape", ] <- c(0, 0, 1)
    }

    ## transformation jacobian
    if ( !(distname.y %in% c("gpd", "GPD")) ) {
      ##transNames <- names(transParRen)
      ##cat("transNames = ", transNames, "\n")
      derTrans  <- matrix(0, nrow = 4,  ncol = 4)
      rownames(derTrans) <- c("lambda", "threshold", transNames)
      colnames(derTrans) <- c("lambda", "threshold", "scale", "shape")
      derTrans["lambda", "lambda"] <- 1.0
      derTrans["threshold", "threshold"] <- 1.0
      derTrans[transNames, c("scale", "shape")] <-
        attr(transParRen, "jacobian")
      derParRen <- derTrans %*% derParRen 
    }

    attr(parRen, "jacobian") <- derParRen

  }

  if (!is.null(vcovGev)) {
    if ( !is.matrix(vcovGev) || (nrow(vcovGev) != 3L) || (ncol(vcovGev) != 3L) ) {
      stop("'vcovGev' must be a matrix with 3 rows and 3 columns")
    }
    if ( !setequal(rownames(vcovGev), gevNames) ) {
      stop("'vcovGev' must have suitable rownames and colnames, see help")
    }
    colnames(vcovGev) <- rownames(vcovGev)
    vcovRen <- derParRen %*% vcovGev[gevNames,gevNames] %*% t(derParRen)
    attr(parRen, "vcov") <- vcovRen
  } 
  
  attr(parRen, "threshold") <- as.numeric(threshold)
  attr(parRen, "distname.y") <- as.character(distname.y)
  
  parRen    
  
}
