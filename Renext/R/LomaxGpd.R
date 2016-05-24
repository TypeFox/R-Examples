##***********************************************************************
## Parameter transfgormation and its jacobian for 
##     
##  Lomax <-> GPD (shape > 0),  maxlo <-> GPD (shape < 0)    
##
## Author: Yves Deville <deville.yves@alpestat.com> 
##
## Note that the lamx and maxlo distribution may accept a 'loc'
## parameter in the future.
##
##***********************************************************************

gpd2lomax <- function(parGpd,
                      jacobian = TRUE,
                      vcovGpd = NULL) {

  gpdNames <- c("scale", "shape")
  lomaxNames <- c("scale", "shape")

  ## checks
  if ( !setequal(names(parGpd), gpdNames) ) {
    stop("'parGpd' must have suitable names, see help")
  }
  if (parGpd["shape"] <= 0){
    stop("\"gpd\" to \"lomax\" conversion possible only when shape > 0")
  }

  ## transform
  parLomax <- c("scale" = as.numeric(parGpd["scale"] / parGpd["shape"]),
                "shape" = as.numeric(1 / parGpd["shape"]))

  if (jacobian) {
    derParLomax <- matrix(0, nrow = 2L, ncol = 2L)
    rownames(derParLomax) <- lomaxNames
    colnames(derParLomax) <- gpdNames
    derParLomax["scale", "scale"] <- 1.0 / parGpd["shape"]
    derParLomax["scale", "shape" ] <- -parGpd["scale"] / (parGpd["shape"]^2)
    derParLomax["shape", "scale"] <- 0.0
    derParLomax["shape", "shape" ] <- -1.0 / (parGpd["shape"]^2)
    attr(parLomax, "jacobian") <- derParLomax
  }
  
  if (!is.null(vcovGpd)) {
    if ( !is.matrix(vcovGpd) || (nrow(vcovGpd) != 2L)
        || (ncol(vcovGpd) != 2L) ) {
      stop("'vcovGpd' must be a matrix with 2 rows and 2 columns")
    }
    if ( !setequal(rownames(vcovGpd), gpdNames) ) {
      stop("'vcovGpd' must have suitable rownames and colnames, see help")
    }
    colnames(vcovGpd) <- rownames(vcovGpd)
    vcovLomax <-
      derParLomax %*% vcovGpd[gpdNames, gpdNames] %*% t(derParLomax)
    attr(parLomax, "vcov") <- vcovLomax
  }
  
  parLomax
  
}

##=======================================================================
## Lomax to gpd
##=======================================================================

lomax2gpd <- function(parLomax,
                      jacobian = TRUE,
                      vcovLomax = NULL) {

  gpdNames <- c("scale", "shape")
  lomaxNames <- c("scale", "shape")
  
  if ( !setequal(names(parLomax), lomaxNames) ) {
    stop("'parLomax' must have suitable names, see help")
  }
  if (parLomax["scale"] <= 0 || parLomax["shape"] <= 0 ){
    stop("bad Lomax parameters on entry")
  }

  parGpd <- c("scale" = as.numeric(parLomax["scale"] / parLomax["shape"]),
              "shape" = as.numeric(1.0 / parLomax["shape"]))
  
  if (jacobian) {
    derParGpd <- matrix(0, nrow = 2L, ncol = 2L)
    rownames(derParGpd) <- gpdNames
    colnames(derParGpd) <- lomaxNames
    derParGpd["scale", "scale"] <- 1.0 / parLomax["shape"]
    derParGpd["scale", "shape" ] <- -parLomax["scale"] / (parLomax["shape"]^2)
    derParGpd["shape", "scale"] <- 0.0
    derParGpd["shape",  "shape" ] <- -1.0 / (parLomax["shape"]^2)
    attr(parGpd, "jacobian") <- derParGpd
  }
  
  if (!is.null(vcovLomax)) {
    if ( !is.matrix(vcovLomax) || (nrow(vcovLomax) != 2L) ||
        (ncol(vcovLomax) != 2L) ) {
      stop("'vcovLomax' must be a matrix with 2 rows and 2 columns")
    }
    if ( !setequal(rownames(vcovLomax), lomaxNames) ) {
      stop("'vcovLomax' must have suitable rownames and colnames, see help")
    }
    colnames(vcovLomax) <- rownames(vcovLomax)
    vcovGpd <-
      derParGpd %*% vcovLomax[lomaxNames, lomaxNames] %*% t(derParGpd)
    attr(parGpd, "vcov") <- vcovGpd
  } 

  parGpd
  
}

##=======================================================================
## gpd to maxlo
##=======================================================================

gpd2maxlo <- function(parGpd,
                      jacobian = TRUE,
                      vcovGpd = NULL) {
  
  gpdNames <- c("scale", "shape")
  maxloNames <- c("scale", "shape")

  ## checks
  if ( !setequal(names(parGpd), gpdNames) ) {
    stop("'parGpd' must have suitable names, see help")
  }
  if (parGpd["shape"] >= 0){
    stop("\"gpd\" to \"maxlo\" conversion possible only when shape > 0")
  }

  ## transform
  parMaxlo <- c("scale" = as.numeric(-parGpd["scale"] / parGpd["shape"]),
                "shape" = as.numeric(-1 / parGpd["shape"]))

  if (jacobian) {
    derParMaxlo <- matrix(0, nrow = 2L, ncol = 2L)
    rownames(derParMaxlo) <- maxloNames
    colnames(derParMaxlo) <- gpdNames
    derParMaxlo["scale", "scale"] <- -1.0 / parGpd["shape"]
    derParMaxlo["scale", "shape" ] <- parGpd["scale"] / (parGpd["shape"]^2)
    derParMaxlo["shape", "scale"] <- 0.0
    derParMaxlo["shape", "shape" ] <- 1.0 / (parGpd["shape"]^2)
    attr(parMaxlo, "jacobian") <- derParMaxlo
  }
  
  if (!is.null(vcovGpd)) {
    if ( !is.matrix(vcovGpd) || (nrow(vcovGpd) != 2L) ||
        (ncol(vcovGpd) != 2L) ) {
      stop("'vcovGpd' must be a matrix with 2 rows and 2 columns")
    }
    if ( !setequal(rownames(vcovGpd), gpdNames) ) {
      stop("'vcovGpd' must have suitable rownames and colnames, see help")
    }
    colnames(vcovGpd) <- rownames(vcovGpd)
    vcovMaxlo <-
      derParMaxlo %*% vcovGpd[gpdNames, gpdNames] %*% t(derParMaxlo)
    attr(parMaxlo, "vcov") <- vcovMaxlo
  }
  
  parMaxlo
  
}

##=======================================================================
## Maxlo to gpd
##=======================================================================

maxlo2gpd <- function(parMaxlo,
                      jacobian = TRUE,
                      vcovMaxlo = NULL) {

  gpdNames <- c("scale", "shape")
  maxloNames <- c("scale", "shape")

  if (parMaxlo["shape"] <= 2.0) {
    if (!is.null(vcovMaxlo)) {
      warning("in 'parMaxlo', 'shape' is <= 2. 'vcovMaxLo' ignored")
      vcovMaxlo <- NULL
    }
  }
  
  if ( !setequal(names(parMaxlo), maxloNames) ) {
    stop("'parMaxlo' must have suitable names, see help")
  }
  if (parMaxlo["scale"] <= 0 || parMaxlo["shape"] <= 0 ){
    stop("bad Maxlo parameters on, entry")
  }

  parGpd <- c("scale" = as.numeric(parMaxlo["scale"] / parMaxlo["shape"]),
              "shape" = as.numeric(-1.0 / parMaxlo["shape"]))
  
  if (jacobian) {
    derParGpd <- matrix(0, nrow = 2L, ncol = 2L)
    rownames(derParGpd) <- gpdNames
    colnames(derParGpd) <- maxloNames
    derParGpd["scale", "scale"] <- 1.0 / parMaxlo["shape"]
    derParGpd["scale", "shape" ] <- -parMaxlo["scale"] / (parMaxlo["shape"]^2)
    derParGpd["shape", "scale"] <- 0.0
    derParGpd["shape",  "shape" ] <- 1.0 / (parMaxlo["shape"]^2)
    attr(parGpd, "jacobian") <- derParGpd
  }
  
  if (!is.null(vcovMaxlo)) {
    if ( !is.matrix(vcovMaxlo) || (nrow(vcovMaxlo) != 2L) ||
        (ncol(vcovMaxlo) != 2L) ) {
      stop("'vcovMaxlo' must be a matrix with 2 rows and 2 columns")
    }
    if ( !setequal(rownames(vcovMaxlo), maxloNames) ) {
      stop("'vcovMaxlo' must have suitable rownames and colnames, see help")
    }
    colnames(vcovMaxlo) <- rownames(vcovMaxlo)
    vcovGpd <-
      derParGpd %*% vcovMaxlo[maxloNames, maxloNames] %*% t(derParGpd)
    attr(parGpd, "vcov") <- vcovGpd
  } 

  parGpd
  
}


