# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## workhorse function for groupwise LARS
#' @import parallel

grouplars <- function(x, y, sMax = NA, assign, robust = FALSE, 
                      centerFun = mean, scaleFun = sd, 
                      regFun = lm.fit, regArgs = list(), 
                      combine = c("min", "euclidean", "mahalanobis"), 
                      winsorize = FALSE, const = 2, prob = 0.95, 
                      fit = TRUE, s = c(0, sMax), crit = c("BIC", "PE"), 
                      splits = foldControl(), cost = rmspe, costArgs = list(), 
                      selectBest = c("hastie", "min"), seFactor = 1, 
                      ncores = 1, cl = NULL, seed = NULL, model = TRUE, 
                      call = NULL) {
  ## initializations
  n <- length(y)
  assignList <- split(seq_len(length(assign)), assign)  # indices per group
  m <- length(assignList)          # number of groups
  p <- sapply(assignList, length)  # number of variables in each group
  adjust <- length(unique(p)) > 1  # adjust for group length?
  if(isTRUE(is.numeric(sMax) && length(sMax) == 2)) {
    # ensure backwards compatibility concerning the number of steps
    s <- c(0, sMax[2])
    sMax <- s[1]
  }
  robust <- isTRUE(robust)
  sMax <- checkSMax(sMax, n, p, groupwise=TRUE)  # number of groups to sequence
  if(robust) {
    # if possible, do not use formula interface
    regControl <- getRegControl(regFun)
    regFun <- regControl$fun
    callRegFun <- getCallFun(regArgs)
    # check how data cleaning weights should be combined
    combine <- match.arg(combine)
    # alternatively use winsorization
    winsorize <- isTRUE(winsorize)
  }
  if(is.na(ncores)) ncores <- detectCores()  # use all available cores
  if(!is.numeric(ncores) || is.infinite(ncores) || ncores < 1) {
    ncores <- 1  # use default value
    warning("invalid value of 'ncores'; using default value")
  } else ncores <- as.integer(ncores)
  
  ## check whether submodels along the sequence should be computed
  ## if yes, check whether the final model should be found via resampling-based 
  ## prediction error estimation
  fit <- isTRUE(fit)
  if(fit) {
    s <- checkSRange(s, sMax=sMax)  # check range of steps along the sequence
    crit <- if(!is.na(s[2]) && s[1] == s[2]) "none" else match.arg(crit)
    if(crit == "PE") {
      # further checks for steps along the sequence
      if(is.na(s[2])) {
        s[2] <- min(sMax, floor(n/2))
        if(s[1] > sMax) s[1] <- sMax
      }
      # set up function call to be passed to perryFit()
      remove <- c("x", "y", "crit", "splits", "cost", "costArgs", 
                  "selectBest", "seFactor", "ncores", "cl", "seed")
      remove <- match(remove, names(call), nomatch=0)
      funCall <- call[-remove]
      funCall$sMax <- sMax
      funCall$s <- s
      # call function perryFit() to perform prediction error estimation
      s <- seq(from=s[1], to=s[2])
      selectBest <- match.arg(selectBest)
      out <- perryFit(funCall, x=x, y=y, splits=splits, 
                      predictArgs=list(s=s, recycle=TRUE), cost=cost, 
                      costArgs=costArgs, envir=parent.frame(2), 
                      ncores=ncores, cl=cl, seed=seed)
      out <- perryReshape(out, selectBest=selectBest, seFactor=seFactor)
      fits(out) <- s
      # fit final model
      funCall$x <- call$x
      funCall$y <- call$y
      funCall$s <- s[out$best]
      funCall$ncores <- call$ncores
      funCall$cl <- cl
      out$finalModel <- eval(funCall, envir=parent.frame(2))
      out$call <- call
      # assign class and return object
      class(out) <- c("perrySeqModel", class(out))
      return(out)
    }
  }
  
  ## prepare the data
  if(!is.null(seed)) set.seed(seed)
  if(fit || (robust && !winsorize)) {
    # check whether parallel computing should be used
    haveCl <- inherits(cl, "cluster")
    haveNcores <- !haveCl && ncores > 1
    useParallel <- haveNcores || haveCl
    # set up multicore or snow cluster if not supplied
    if(haveNcores) {
      if(.Platform$OS.type == "windows") {
        cl <- makePSOCKcluster(rep.int("localhost", ncores))
      } else cl <- makeForkCluster(ncores)
      on.exit(stopCluster(cl))
    }
    if(useParallel) {
      # set seed of the random number stream
      if(!is.null(seed)) clusterSetRNGStream(cl, iseed=seed)
      else if(haveNcores) clusterSetRNGStream(cl)
    }
  }
  if(robust) {
    # use fallback mode (standardization with mean/SD) for dummies 
    z <- robStandardize(y, centerFun, scaleFun)
    xs <- robStandardize(x, centerFun, scaleFun, fallback=TRUE)
    muY <- attr(z, "center")
    sigmaY <- attr(z, "scale")
    muX <- attr(xs, "center")
    sigmaX <- attr(xs, "scale")
    if(winsorize) {
      if(is.null(const)) const <- 2
      # obtain data cleaning weights from winsorization
      w <- winsorize(cbind(z, xs), standardized=TRUE, const=const, 
                     prob=prob, return="weights")
    } else {
      # clean data in a limited sense: there may still be correlation 
      # outliers between the groups, but these should not be a problem
      # compute weights from robust regression for each group
      # intercept column needs to be used in robust short regressions
      if(combine == "min") {
        # define function to compute weights for each predictor group
        getWeights <- function(i, x, y) {
          x <- x[, i, drop=FALSE]
          if(regControl$useFormula) {
            fit <- callRegFun(y ~ x, fun=regFun, args=regArgs)
          } else {
            x <- cbind(rep.int(1, n), x)
            fit <- callRegFun(x, y, fun=regFun, args=regArgs)
          }
          weights(fit)
        }
        # compute weights for each predictor group
        if(useParallel) {
          w <- parSapply(cl, assignList, getWeights, xs, z)
        } else w <- sapply(assignList, getWeights, xs, z)
        w <- sqrt(apply(w, 1, min))  # square root of smallest weights
        # observations can have zero weight, in which case the number 
        # of observations needs to be adjusted
        n <- length(which(w > 0))
      } else {
        # define function to compute scaled residuals for each 
        # predictor group
        getResiduals <- function(i, x, y) {
          x <- x[, i, drop=FALSE]
          if(regControl$useFormula) {
            fit <- callRegFun(y ~ x, fun=regFun, args=regArgs)
          } else {
            x <- cbind(rep.int(1, n), x)
            fit <- callRegFun(x, y, fun=regFun, args=regArgs)
          }
          residuals(fit) / getScale(fit)
        }
        # compute scaled residuals for each predictor group
        if(useParallel) {
          residuals <- parSapply(cl, assignList, getResiduals, xs, z)
        } else residuals <- sapply(assignList, getResiduals, xs, z)
        # obtain weights from scaled residuals
        if(combine == "euclidean") {
          # assume diagonal structure of the residual correlation matrix
          # and compute weights based on resulting mahalanobis distances
          d <- qchisq(prob, df=m)  # quantile of the chi-squared distribution
          w <- pmin(sqrt(d/rowSums(residuals^2)), 1)
        } else {
          # obtain weights from multivariate winsorization
          w <- winsorize(residuals, standardized=TRUE, const=const, 
                         prob=prob, return="weights")
        }
      }
    }
    z <- standardize(w*z)  # standardize cleaned response
    xs <- standardize(w*xs)
    # center and scale of response
    muY <- muY + attr(z, "center")
    sigmaY <- sigmaY * attr(z, "scale")
    # center and scale of candidate predictor variables
    muX <- muX + attr(xs, "center")
    sigmaX <- sigmaX * attr(xs, "scale")
  } else {
    z <- standardize(y)   # standardize response
    xs <- standardize(x)
    # center and scale of response
    muY <- attr(z, "center")
    sigmaY <- attr(z, "scale")
    # center and scale of candidate predictor variables
    muX <- attr(xs, "center")
    sigmaX <- attr(xs, "scale")
  }
  
#   ## find first ranked predictor group
#   zHat <- sapply(assignList, 
#                  function(i, x, y) {
#                    x <- x[, i, drop=FALSE]
#                    model <- lm.fit(x, y)
#                    fitted(model)
#                  }, xs, z)
#   corZ <- unname(apply(zHat, 2, sd))
#   # if necessary, compute the denominators for adjustment w.r.t. the number
#   # of variables in the group, and adjust the correlations of the fitted
#   # values with the response
#   if(adjust) {
#     adjustment <- sqrt(p)
#     corZ <- corZ / adjustment
#   }
#   # first active group
#   active <- which.max(corZ)
#   r <- c(corZ[active], rep.int(NA, sMax[1]-1))
#   # not yet sequenced groups
#   inactive <- seq_len(m)[-active]
#   corZ <- corZ[-active]
#   
#   ## update active set
#   R <- diag(1, sMax[1])
#   a <- 1
#   if(adjust) a <- a / adjustment[active]
#   w <- 1
#   sigma <- 1
#   # start iterative computations
#   for(k in seq_len(sMax[1]-1)) {
#     # standardize fitted values for k-th group
#     zHat[, active[k]] <- standardize(zHat[, active[k]])
#     # compute the equiangular vector
#     if(k == 1) u <- zHat[, active] # equiangular vector equals first direction
#     else {
#       # compute correlations between fitted values for active groups
#       R[k, seq_len(k-1)] <- R[seq_len(k-1), k] <- 
#         apply(zHat[, active[seq_len(k-1)], drop=FALSE], 2, cor, zHat[, active[k]])
#       # other computations according to algorithm
#       invR <- solve(R[seq_len(k), seq_len(k), drop=FALSE])
#       q <- if(adjust) adjustment[active] else rep.int(1, k)
#       # compute the correlation of the fitted values from the active
#       # predictor groups with the equiangular vector (adjustment for
#       # unequal group size is considered via vector 'q')
#       a <- c(t(q) %*% invR %*% q)^(-1/2)
#       # compute the equiangular vector
#       w <- a * c(invR %*% q)
#       u <- c(zHat[, active, drop=FALSE] %*% w)  # equiangular vector
#     }
#     # compute the fitted values of the equiangular vector for each
#     # inactive predictor group, as well as the correlations involving the
#     # inactive predictor groups and the equiangular vector
#     uHat <- sapply(assignList[inactive], 
#                    function(i, x, u) {
#                      x <- x[, i, drop=FALSE]
#                      model <- lm.fit(x, u)
#                      fitted(model)
#                    }, xs, u)
#     corU <- apply(zHat[, inactive, drop=FALSE], 2, cor, u)
#     tau <- apply(uHat, 2, sd)
#     # adjustment for unequal group size if necessary
#     if(adjust) {
#       tmp <- adjustment[inactive]
#       corU <- corU / tmp
#       tau <- tau / tmp
#     }
#     # compute the step size by solving the quadratic equation
#     gammas <- findStepSizes(r[k], a, corZ, corU, tau)
#     whichMin <- which.min(gammas)
#     gamma <- gammas[whichMin]
#     # the following computations are not necessary in the last iteration
#     if(k < sMax[1]-1) {
#       # update the scale of the current response
#       sigma <- sqrt(1 - 2*gamma*r[k]/a + gamma^2)
#       # update the fitted values from the new or not yet sequenced
#       # predictor groups
#       zHat[, inactive] <- (zHat[, inactive] - gamma * uHat) / sigma
#       r[k+1] <- (r[k] - gamma * a) / sigma
#       corZ <- corZ[-whichMin]
#       corU <- corU[-whichMin]
#       tau <- tau[-whichMin]
#       corZ = sqrt(corZ^2 - 2*gamma*corU*corZ + gamma^2 * tau^2) / sigma
#     }
#     # update active set and not yet sequenced groups
#     active <- c(active, inactive[whichMin])  # update active set
#     inactive <- inactive[-whichMin]          # update not yet sequenced groups
#   }
  
  ## call C++ function
  active <- .Call("R_fastGrplars", R_x=xs, R_y=z, R_sMax=sMax, 
                  R_assign=assignList, R_ncores=ncores, 
                  PACKAGE="robustHD")
  
  ## choose optimal model according to specified criterion
  if(fit) {
    # add ones to matrix of predictors to account for intercept
    x <- addIntercept(x)
    # call function to fit models along the sequence
    out <- seqModel(x, y, active=active, sMin=s[1], sMax=s[2], 
                    assign=assignList, robust=robust, regFun=regFun, 
                    useFormula=regControl$useFormula, regArgs=regArgs, 
                    crit=crit, cl=cl)
    # add center and scale estimates
    out[c("muX", "sigmaX", "muY", "sigmaY")] <- list(muX, sigmaX, muY, sigmaY)
    # add model data to result if requested
    if(isTRUE(model)) out[c("x", "y")] <- list(x=x, y=y)
    out$assign <- assign  # this is helpful for some methods
    out$robust <- robust
    if(robust) out$w <- w  # add data cleaning weights
    out$call <- call       # add call to return object
    class(out) <- c("grplars", class(out))
    out
  } else active
}


# # find possible step sizes for groupwise LARS by solving quadratic equation
# findStepSizes <- function(r, a, corY, corU, tau) {
#   mapply(function(corY, corU, tau, r, a) {
#     # quadratic equation to be solved
#     comp <- c(r^2 - corY^2, 2 * (corU*corY - a*r), a^2 - tau^2)
#     # solution of quadratic equation
#     gamma <- Re(polyroot(comp))
#     solution <- min(gamma[gamma >= 0])
#     solution
#   }, corY, corU, tau, MoreArgs=list(r=r, a=a))
# }
