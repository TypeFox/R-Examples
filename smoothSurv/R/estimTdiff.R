######################################################################################################
#### AUTHOR:    Arnost Komarek                                                                    ####
####            25/02/2004                                                                        ####
###             03/05/2004                                                                        ####
###             14/01/2010:  debugging for the case that there is only one covariate in a model   ####
###             24/04/2014:  correction to make it work with covariates for the scale parameter   ####
###             28/06/2014:  additional adjustment to keep also covariate labels in the output    ####
####                                                                                              ####
#### FILE:      estimTdiff.R                                                                      ####
####                                                                                              ####
#### FUNCTIONS: estimTdiff                                                                        ####
####            estimTdiff.smoothSurvReg                                                          ####
####                                                                                              ####
######################################################################################################

estimTdiff <- function(x, ...)
{
  UseMethod("estimTdiff")
}


estimTdiff.smoothSurvReg <- function(x, cov1, cov2, logscale.cov1, logscale.cov2, time0 = 0, conf.level=0.95, ...)
{
  if (x$fail >= 99) {
    cat("No estimate of T difference, smoothSurvReg failed.\n")
    return(invisible(x))
  }
  is.intercept <- x$estimated["(Intercept)"]
  common.logscale <- x$estimated["common.logscale"]
  est.scale <- x$estimated["Scale"]
  allregrname <- row.names(x$regres)

## Rownames of cov1 and cov2
## added on 20140628
## =================================
  if (!missing(cov1)){
    if (is.null(rownames(cov1))) rncov1 <- NULL else rncov1 <- rownames(cov1)
  }else{
    rncov1 <- NULL
  }    
  if (!missing(cov2)){
    if (is.null(rownames(cov2))) rncov2 <- NULL else rncov2 <- rownames(cov2)
  }else{
    rncov2 <- NULL
  }    
  
  if (!missing(logscale.cov1)){
    if (is.null(rownames(logscale.cov1))) rnlscov1 <- NULL else rnlscov1 <- rownames(logscale.cov1)
  }else{
    rnlscov1 <- NULL
  }    
  if (!missing(logscale.cov2)){
    if (is.null(rownames(logscale.cov2))) rnlscov2 <- NULL else rnlscov2 <- rownames(logscale.cov2)
  }else{
    rnlscov2 <- NULL
  }    

  
## INTERCEPT AND SCALE (if it is common)
## =====================================
  mu0 <- ifelse(is.intercept, x$regres["(Intercept)", "Value"], 0)
  if (common.logscale){
    if (est.scale) s0 <- x$regres["Scale", "Value"]
    else           s0 <- x$init.regres["Scale", "Value"]
  }


## COVARIATES FOR REGRESSION
## =========================
  nx <- ncol(x$x)
  ncov <- ifelse(is.intercept, nx - 1, nx)

  ## Manipulate with covariate values from the user
  if (missing(cov1) && missing(cov2) && ncov > 0){
    cov1 <- matrix(rep(0, ncov), nrow = 1)
    cov2 <- cov1
  }
  else
    if((missing(cov1) || missing(cov2)) && ncov > 0) stop("cov1 and cov2 must be either both missing or both present.")

  if (ncov == 0){                      ## only intercept in the model
    cov1 <- NULL
    cov2 <- NULL
  }    
  if (ncov == 1){
    cov1 <- matrix(cov1, ncol = 1)
    cov2 <- matrix(cov2, ncol = 1)
    if (length(cov1) != length(cov2)) stop("cov1 and cov2 must be of same length.")
  }
  if (ncov > 1){
    if (length(cov1) != length(cov2)) stop("cov1 and cov2 must be of same length.")
  }    
  
  ## Different covariates combinations
  row.cov <- ifelse(is.null(dim(cov1)), 1, dim(cov1)[1])
  col.cov <- ifelse(is.null(dim(cov1)),
                    ifelse(is.null(cov1), 0, length(cov1)),
                    dim(cov1)[2])


## COVARIATES FOR LOG-SCALE
## ========================
  nz <- ncol(x$z)
  if (!common.logscale){
    is.intercept.inscale <- (allregrname[nx+1] == "LScale.(Intercept)")
    ncovz <- ifelse(is.intercept.inscale, nz - 1, nz)

    ## logscale: Manipulate with covariate values from the user
    if (missing(logscale.cov1) && missing(logscale.cov2) && ncovz > 0){
      logscale.cov1 <- matrix(rep(0, ncovz), nrow = 1)
      logscale.cov2 <- logscale.cov1
    }
    else
      if((missing(logscale.cov1) || missing(logscale.cov2)) && ncovz > 0) 
        stop("logscale.cov1 and logscale.cov2 must be either both missing or both present.")

    if (ncovz == 0){                      ## only intercept in the model for log-scale (this should never happen if !common.logscale)
      logscale.cov1 <- NULL
      logscale.cov2 <- NULL
    }    
    if (ncovz == 1){
      logscale.cov1 <- matrix(logscale.cov1, ncol = 1)
      logscale.cov2 <- matrix(logscale.cov2, ncol = 1)
      if (length(logscale.cov1) != length(logscale.cov2)) 
        stop("logscale.cov1 and logscale.cov2 must be of same length.")
    }
    if (ncovz > 1){
      if (length(logscale.cov1) != length(logscale.cov2)) 
        stop("logscale.cov1 and logscale.cov2 must be of same length.")
    }    
  
    ## logscale: Different covariates combinations
    logscale.row.cov <- ifelse(is.null(dim(logscale.cov1)), 1, dim(logscale.cov1)[1])
    logscale.col.cov <- ifelse(is.null(dim(logscale.cov1)),
                               ifelse(is.null(logscale.cov1), 0, length(logscale.cov1)),
                               dim(logscale.cov1)[2])
  }
  else{
    ncovz <- 0
    logscale.row.cov <- row.cov
    logscale.col.cov <- 1
  }    


## LINEAR PREDICTORS
## =================
  beta <- x$regres[1:nx, "Value"]
  if (col.cov != ncov) stop("Incorrect cov1 or cov2 parameter ")
  if (ncov > 0){
    if (is.intercept){
      beta <- matrix(beta[2:nx], nrow = ncov, ncol = 1)
      regrname <- allregrname[2:nx]
    }
    else{
      beta <- matrix(beta[1:nx], nrow = ncov, ncol = 1)
      regrname <- allregrname[1:nx]
    }
    cov1 <- matrix(cov1, nrow = row.cov, ncol = col.cov)
    cov2 <- matrix(cov2, nrow = row.cov, ncol = col.cov)    
    eta1 <- mu0 + as.numeric(cov1 %*% beta)
    eta2 <- mu0 + as.numeric(cov2 %*% beta)
    if (is.intercept){
      cov1 <- cbind(rep(1, row.cov), cov1)               ## add intercept to the covariates
      cov2 <- cbind(rep(1, row.cov), cov2)               ## to be used in further computation
    }
  }
  else{    ## only intercept in the model
     regrname <- character(0)
     eta1 <- rep(mu0, row.cov)
     eta2 <- rep(mu0, row.cov)
     cov1 <- matrix(1, nrow = row.cov, ncol = 1)
     cov2 <- cov1
  }

  
## LINEAR PREDICTORS FOR LOG-SCALE, AND COMPUTATION OF A SCALE
## ===========================================================
  if (!common.logscale){
    pars.scale <- x$regres[(nx+1):(nx+nz), "Value"]
    if (logscale.col.cov != ncovz) stop("Incorrect logscale.cov1 or logscale.cov2 parameter ")
    if (row.cov != logscale.row.cov) stop("Different number of covariate combinations for regression and log-scale ")

    if (ncovz > 0){
      if (is.intercept.inscale){
        sint <- pars.scale[1]
        pars.scale <- matrix(pars.scale[2:nz], nrow = ncovz, ncol = 1)
        scalename <- allregrname[(nx+2):(nx+nz)]
      }
      else{
        sint <- 0
        pars.scale <- matrix(pars.scale[1:nz], nrow = ncovz, ncol = 1)
        scalename <- allregrname[(nx+1):(nx+nz)]
      }
      logscale.cov1 <- matrix(logscale.cov1, nrow = logscale.row.cov, ncol = logscale.col.cov)
      logscale.cov2 <- matrix(logscale.cov2, nrow = logscale.row.cov, ncol = logscale.col.cov)
      logscale1 <- sint + as.numeric(logscale.cov1 %*% pars.scale)
      logscale2 <- sint + as.numeric(logscale.cov2 %*% pars.scale)
      if (is.intercept.inscale){
        logscale.cov1 <- cbind(rep(1, row.cov), logscale.cov1)      ## add intercept to the covariates
        logscale.cov2 <- cbind(rep(1, row.cov), logscale.cov2)      ## to be used in further computation
      }
    }
    else{    ## this should never happen if !common.logscale
       sint <- pars.scale[1]
       scalename <- character(0)
       logscale1 <- rep(sint, logscale.row.cov)
       logscale2 <- rep(sint, logscale.row.cov)
       logscale.cov1 <- matrix(1, nrow = logscale.row.cov, ncol = 1)
       logscale.cov2 <- logscale.cov1
    }
    s0.1 <- exp(logscale1)
    s0.2 <- exp(logscale2)
  }
  else{
    scalename <- character(0)    
    s0.1 <- rep(s0, row.cov)
    s0.2 <- s0.1
    logscale.cov1 <- matrix(1, nrow = logscale.row.cov, ncol = 1)
    logscale.cov2 <- logscale.cov1    
  }


## ESTIMATION AND STANDARD ERRORS
## ==============================
  ## Compute M = E(s0*epsilon) (for each covariate combination)
  knots <- as.numeric(x$spline$Knot)
  ccoef <- as.numeric(x$spline[["c coef."]])
  sigmasq0 <- (x$spline["SD basis"])^2  
  nknots <- length(knots)

  s0sq.1 <- s0.1^2
  s0sq.2 <- s0.2^2

  ccoefrep <- matrix(rep(ccoef, row.cov), ncol = nknots, byrow = TRUE)
  knotsrep <- matrix(rep(knots, row.cov), ncol = nknots, byrow = TRUE)
  sigmasq0rep <- matrix(rep(sigmasq0, row.cov), ncol = nknots, byrow = TRUE)

  s0.knots.1 <- matrix(s0.1 %x% knots, ncol = nknots, byrow = TRUE)
  s0sq.sigmasq0.1 <- matrix(s0sq.1 %x% sigmasq0, ncol = nknots, byrow = TRUE)
  expm.1 <- exp(s0.knots.1 + 0.5*s0sq.sigmasq0.1)
  c.expm.1 <- ccoefrep * expm.1
  M.1 <- apply(c.expm.1, 1, sum)
  
  s0.knots.2 <- matrix(s0.2 %x% knots, ncol = nknots, byrow = TRUE)
  s0sq.sigmasq0.2 <- matrix(s0sq.2 %x% sigmasq0, ncol = nknots, byrow = TRUE)
  expm.2 <- exp(s0.knots.2 + 0.5*s0sq.sigmasq0.2)
  c.expm.2 <- ccoefrep * expm.2
  M.2 <- apply(c.expm.2, 1, sum)

  ## Estimates of E(T) and a difference
  expeta1 <- exp(eta1)
  expeta2 <- exp(eta2)
  ET1 <- M.1 * expeta1
  ET2 <- M.2 * expeta2
  diffT <- ET1 - ET2

  ## Derivatives of estimates w.r.t. regression
  ## (at each column is a derivative of a quantity for one specific covar. value)
  dET1 <- t(cov1) * matrix(ET1, nrow = ncov + 1, ncol = row.cov, byrow = TRUE)
  dET2 <- t(cov2) * matrix(ET2, nrow = ncov + 1, ncol = row.cov, byrow = TRUE)
  ddiffT <- dET1 - dET2

  ## Derivatives of E(s0*epsilon)
  if (est.scale){
    c.mu.expm.1 <- knotsrep * expm.1
    c.sigmasq0.expm.1 <- sigmasq0rep * expm.1
    s0.1rep <- matrix(rep(s0.1, nknots), ncol = nknots)    
#    dMdgamma.1 <- apply(s0.1rep*c.mu.expm.1 + s0.1rep^2*c.sigmasq0.expm.1, 1, sum)
    dMdgamma.1 <- apply(s0.1rep^2*c.sigmasq0.expm.1, 1, sum)    
    dET1dgamma <- dMdgamma.1 * exp(eta1)    
    
    c.mu.expm.2 <- knotsrep * expm.2
    c.sigmasq0.expm.2 <- sigmasq0rep * expm.2
    s0.2rep <- matrix(rep(s0.2, nknots), ncol = nknots)    
#    dMdgamma.2 <- apply(s0.2rep*c.mu.expm.2 + s0.2rep^2*c.sigmasq0.expm.2, 1, sum)
    dMdgamma.2 <- apply(s0.2rep^2*c.sigmasq0.expm.2, 1, sum)    
    dET2dgamma <- dMdgamma.2 * exp(eta2)        

    if (!common.logscale){
      dET1ds <- t(logscale.cov1) * matrix(dET1dgamma, nrow = ncol(logscale.cov1), ncol = row.cov, byrow = TRUE)
      dET2ds <- t(logscale.cov2) * matrix(dET2dgamma, nrow = ncol(logscale.cov2), ncol = row.cov, byrow = TRUE)
      dET1 <- rbind(dET1, dET1ds)
      dET2 <- rbind(dET2, dET2ds)
      ddiffT <- rbind(ddiffT, dET1ds - dET2ds)      
    }      
    else{
      dET1 <- rbind(dET1, dET1dgamma)
      dET2 <- rbind(dET2, dET2dgamma)
      ddiffT <- rbind(ddiffT, dET1dgamma - dET2dgamma)      
    }      
  }    

  if (x$estimated["ccoef"]){
    dMdD.1 <- matrix(0, nrow = nknots - 3, ncol = row.cov)
    dMdD.2 <- matrix(0, nrow = nknots - 3, ncol = row.cov)    
    for (i in 1:row.cov){
      dMdD.1[, i] <- matrix(apply(matrix(expm.1[i, ], ncol = nknots, nrow = nknots - 3, byrow = TRUE) * x$dCdD, 1, sum), ncol = 1)
      dMdD.2[, i] <- matrix(apply(matrix(expm.2[i, ], ncol = nknots, nrow = nknots - 3, byrow = TRUE) * x$dCdD, 1, sum), ncol = 1)
    }
    dET1dD <- dMdD.1 * matrix(exp(eta1), nrow = nknots - 3, ncol = row.cov, byrow = TRUE)
    dET2dD <- dMdD.2 * matrix(exp(eta2), nrow = nknots - 3, ncol = row.cov, byrow = TRUE)
    dET1 <- rbind(dET1, dET1dD)
    dET2 <- rbind(dET2, dET2dD)
    ddiffT <- rbind(ddiffT, dET1dD - dET2dD)
  }

  ## Variance estimates for ET1, ET2, ddiffT
  varET1 <- diag(t(dET1) %*% x$var %*% dET1)
  varET2 <- diag(t(dET2) %*% x$var %*% dET2)
  vardiffT <- diag(t(ddiffT) %*% x$var %*% ddiffT)

  sdET1 <- sqrt(varET1)
  sdET2 <- sqrt(varET2)
  sddiffT <- sqrt(vardiffT)

  alpha <- 1 - conf.level
  if (alpha <= 0 | alpha >= 1) stop("incorrect confidence level supplied")
  u <- qnorm(1 - alpha/2)
  ET1.lower <- ET1 + time0 - sdET1*u
  ET1.lower[ET1.lower < 0] <- 0
  ET1.upper <- ET1 + time0 + sdET1*u  
  #
  ET2.lower <- ET2 + time0 - sdET2*u
  ET2.lower[ET2.lower < 0] <- 0
  ET2.upper <- ET2 + time0 + sdET2*u  
  #
  diffT.lower <- diffT - sddiffT*u
  diffT.upper <- diffT + sddiffT*u
    
  obj <- data.frame(ET1 = ET1 + time0, sd.ET1 = sdET1, ET1.lower= ET1.lower, ET1.upper= ET1.upper, 
                    ET2 = ET2 + time0, sd.ET2 = sdET2, ET2.lower= ET2.lower, ET2.upper= ET2.upper, 
                    diffT = diffT, sd.diffT = sddiffT, diffT.lower=diffT.lower, diffT.upper=diffT.upper)
  rownames(obj) <- 1:nrow(obj)

  if (ncov > 0){                 ### 14/01/2010:  changed from 'if (ncov > 1)' which was not correct 
    cov1 <- cov1[, -1]                   ## Remove intercept from the covariates
    cov2 <- cov2[, -1]
    if (ncov == 1){                      ### added on 14/01//2010
      cov1 <- matrix(cov1, ncol=1)       ### added on 14/01//2010
      cov2 <- matrix(cov2, ncol=1)       ### added on 14/01//2010
    }                                    ### added on 14/01//2010
    colnames(cov1) <- regrname
    colnames(cov2) <- regrname
    if (is.null(rncov1)) rownames(cov1) <- paste("Value ", 1:row.cov, sep = "") else rownames(cov1) <- rncov1
    if (is.null(rncov2)) rownames(cov2) <- paste("Value ", 1:row.cov, sep = "") else rownames(cov2) <- rncov2
  }
  else{
    cov1 <- NULL
    cov2 <- NULL
  }    

  if (ncovz > 1){
    if (is.intercept.inscale){
      logscale.cov1 <- logscale.cov1[, -1]                   ## Remove intercept from the covariates
      logscale.cov2 <- logscale.cov2[, -1]
    }      
    colnames(logscale.cov1) <- scalename
    colnames(logscale.cov2) <- scalename
    if (is.null(rnlscov1)) rownames(logscale.cov1) <- paste("Value ", 1:row.cov, sep = "") else rownames(logscale.cov1) <- rnlscov1
    if (is.null(rnlscov2)) rownames(logscale.cov2) <- paste("Value ", 1:row.cov, sep = "") else rownames(logscale.cov2) <- rnlscov2
  }    

  attr(obj, "conf.level") <- conf.level
  attr(obj, "cov1") <- cov1
  attr(obj, "cov2") <- cov2
  if (!common.logscale){
    attr(obj, "logscale.cov1") <- logscale.cov1
    attr(obj, "logscale.cov2") <- logscale.cov2    
  }
  else{
    attr(obj, "logscale.cov1") <- NULL
    attr(obj, "logscale.cov2") <- NULL
  }    
  class(obj) <- "estimTdiff"
  
  return(obj)        
}  
