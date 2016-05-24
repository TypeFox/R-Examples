################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Auxiliary functions for twinSIR()
### and to compute one-sided AIC by simulation (in twinSIR_methods.R)
###
### Copyright (C) 2009-2014 Sebastian Meyer, contributions by Michael Hoehle
### $Revision: 991 $
### $Date: 2014-09-01 23:13:26 +0200 (Mon, 01. Sep 2014) $
################################################################################


################################################################################
# The cox function is used in model formulae to indicate/capture the variables
# which go into the cox part/endemic component of the model.
# Also, with this "cox variables" it is possible to build up interactions
# as usual: cox(var1):cox(var2)... (as if cox(...) was a normal variable)
################################################################################

cox <- function (x)
{
    x
}


################################################################################
# read.design extracts the two parts X and Z of the design matrix.
# Z contains the endemic part (consisting of the cox(.) terms),
# X contains the epidemic part (the rest).
# The automatic intercept variable is excluded from these matrices!
#
# ARGS:
#  m - a model.frame
#  Terms - terms for this model.frame (used to extract the model.matrix from m)
# RETURNS:
#  list of matrices X and Z.
#  If there is no variable in one part of the model the corresponding matrix has
#  0 columns, e.g. ncol(Z) = 0, if there is no endemic (Cox) part.
# NOTE:
# This function is inspired from the timereg package by T. Scheike (available
# under GPL2). See http://staff.pubhealth.ku.dk/~ts/timereg.html for details.
# The function has been improved/modified to fit our purposes.
################################################################################

read.design <- function (m, Terms)
{
    attr(Terms, "intercept") <- 1   # we will remove the intercept later on
    # we need this to ensure that we have a reference category
    # in case of factors (correct contrasts)
    XZ <- model.matrix(Terms, m)
    Zterms <- grep("cox\\([^)]+\\)", colnames(XZ), ignore.case = FALSE,
        perl = FALSE, value = FALSE, fixed = FALSE, useBytes = FALSE,
        invert = FALSE)
    # timereg 1.0-9 way: pattern="^cox[(][A-z0-9._]*[)]" with perl=TRUE
    
    X <- XZ[, -c(1L, Zterms), drop = FALSE]
    Z <- XZ[, Zterms, drop = FALSE]
    
    ud <- list(X = X, Z = Z)
    return(ud)
}

## Alternative way to do the same thing as read.design.
## This approach is similar to that of coxph, but most often some milliseconds
## slower.
# read.design <- function (m, Terms)
# {
#     attr(Terms, "intercept") <- 1   # we will remove the intercept later on
#     # we need this to ensure that we have a reference category
#     # in case of factors (right contrasts)
#     nCoxTerms <- length(attr(Terms, "specials")[["cox"]])
#     if (nCoxTerms > 0) {
#         dropX <- untangle.specials(Terms, "cox", order=1:3)$terms
#     }
#     if (length(dropX) > 0) {
#         X <- model.matrix(Terms[-dropX], m)  # by subscripting a Terms object,
#         Z <- model.matrix(Terms[dropX], m)   # one always gets an intercept term
#         Z <- Z[, -1, drop = FALSE]
#     } else {
#         X <- model.matrix(Terms, m)
#         Z <- X[, NULL, drop = FALSE]
#     }
#     X <- X[, -1, drop = FALSE]
#     
#     ud <- list(X = X, Z = Z)
#     return(ud)
# }


################################################################################
# Little helper function which returns either summary(object) or simply object,
# if it is already a summary. The function also verifies the 'class'.
################################################################################

getSummary <- function (object, class)
{
    summaryClass <- paste("summary", class, sep=".")
    if (inherits(object, class)) {
        summary(object)
    } else if (inherits(object, summaryClass)) {
        object
    } else {
        stop("'object' must inherit from class \"", summaryClass, "\" or \"",
             class, "\"")
    }
}


################################################################################
############################## OSAIC function ##################################
################################################################################
# Two functions:
# Ztilde.chibarsq <- function(Z,p,Winv,tR,s=1)
# w.chibarsq.sim <- function(p, W, N=1e4) 
#
# Both functions are only used internally, no need for documentation
# they are used in function .OSAICpenalty  (twinSIR_methods.R)
################################################################################

##########################################################################
# This function computes Ztilde
# for one Z as specified in Simulation 3, Silvapulle & Sen (2005), p. 79.
# See also p. 37 for the quadprog link.
#
# Params:
#  Z - px1 matrix or vector with specific Z value
#  p - dimension of the problem, where theta is restricted to R^{+p}
#  Winv - inverse of covariance matrix of Z
#  tR - transpose of constraint matrix R\theta \geq 0. In all cases equal to
#      diag(p), but to save time we deliver it to the function every time
#  s - rescale objective function (division by s)
#
# Returns:
#  Ztilde, the point at which (Z-\theta)' W^{-1} (Z-\theta) is the
#  minimum over \theta \geq 0.
##########################################################################

Ztilde.chibarsq <- function(Z,p,Winv,tR,s=1)  {
  #The solve.QP function minimizes
  #-d^T b + 1/2 b^T D b subject to the constraints A^T b >= b_0.
  #Thus using p. 37 we have d = t(Winv) %*% Z.
  d <- crossprod(Winv, Z)
  #Note: Winv and d can become quiet large (or small), but since the solution is
  #invariant to the scaling of the function being minimized, we can equivalently
  #call solve.QP using D/s and d/s (e.g., s=mean(D)) to avoid the error
  #"constraints are inconsistent, no solution!"
  theta <- quadprog::solve.QP(Dmat = Winv/s, dvec = d/s, Amat = tR,
                              bvec = rep.int(0,p), meq = 0)$solution
  return(sum(theta > 0))
}

######################################################################
# Compute OSAIC by simulation weights as described in Silvapulle & Sen
# (2005), Simulation 3, p.79.
#
# Params:
#  p - dimension of the problem, theta is constrained to R^{+p}
#  W - covariance matrix of the chibarsq distribution
#  N - number of simulations to use
#
# Returns:
#  vector of length p+1 containing the weights w_i, i=0, \ldots, p,
#  computed by Monte Carlo simulation
######################################################################

w.chibarsq.sim <- function(p, W, N=1e4) {
  #Draw Z's from multivariate normal distribution with covariance
  #matrix W
  Z <- mvrnorm(N, rep.int(0,p), W)
  if (is.vector(Z)) Z <- t(Z)  # case N==1
  #inverse of W
  Winv <- solve(W)
  #For each simulation calculate Ztilde
  sims <- apply(X=Z, MARGIN=1, FUN=Ztilde.chibarsq,
                p=p, Winv=Winv, tR=diag(p), s=mean(Winv))
  w <- table(factor(sims, levels=0:p)) / N
  return(w)
}



################################################################################
# The helper 'getModel.simEpidata' extracts the model of an object of class
# "simEpidata" similar to the function 'twinSIR' with model = TRUE,
# i.e. a list with components survs, X, Z and weights, where atRiskY == 1.
# The log-baseline h0 is evaluated at start times of intervals only.
# This function is used in function 'intensityPlot'.
################################################################################

getModel.simEpidata <- function (object, ...)
{
    class(object) <- "data.frame"   # avoid use of [.epidata (not necessary here)
    config <- attr(object, "config")
    alpha <- config$alpha
    beta <- config$beta
    atRiskY1 <- object$atRiskY == 1
    simepi1 <- object[atRiskY1,]
    survs <- simepi1[c("id", "start", "stop", "event")]
    attr(survs, "eventTimes") <- attr(object, "eventTimes")
    attr(survs, "timeRange") <- attr(object, "timeRange")
    X <- as.matrix(simepi1[tail(1:ncol(simepi1), length(alpha))])
    logbaseline <- sapply(survs$start, FUN = config$h0, simplify = TRUE)
    Terms <- attr(object, "terms")
    Z <- read.design(model.frame(Terms, simepi1), Terms)$Z
    Z <- cbind("cox(logbaseline)" = logbaseline, Z)
    model <- list(survs = survs, X = X, Z = Z, weights = rep.int(1,nrow(survs)))
    return(model)
}


### Similar auxiliary method extracting the model component
### of a fitted 'twinSIR'

getModel.twinSIR <- function (object, ...)
{
    if (is.null(model <- object[["model"]])) {
        stop("'", deparse(substitute(object)), "' does not contain the 'model' ",
             "component (use 'model = TRUE' when calling 'twinSIR')")
    }
    return(model)
}
