# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## function to compute submodels along a given sequence of predictors
# x ........ predictor matrix including column to account for intercept
# y ........ response
# active ... sequence of predictors
seqModel <- function(x, y, active, sMin = 0, sMax = NA, assign = NULL, 
                     robust = TRUE, regFun = .lmrob.fit, useFormula = FALSE, 
                     regArgs = list(), crit = "BIC", cl = NULL) {
  # initializations
  n <- length(y)
  haveAssign <- !is.null(assign)
  if(haveAssign && !is.list(assign)) {
    # list of column indices for each predictor group
    assign <- split(seq_len(length(assign)), assign)
  }
  if(robust) callRegFun <- getCallFun(regArgs)
  # prepare the variable sequence and the degrees of freedom of the models
  if(haveAssign) {
    # the default is to fit models as long as there are twice as many 
    # observations as predictors
    if(is.na(sMax)) {
      dfMax <- floor(n/2) + 1
      sMax <- dfMax - 1
    } else dfMax <- n
    if(sMax > length(active)) sMax <- length(active)
    s <- 0:sMax
    # compute degrees of freedom of the submodels along sequence
    firstActive <- active[seq_len(sMax)]
    p <- sapply(assign[firstActive], length)  # number of variables per group
    df <- cumsum(c(1, unname(p)))             # degrees of freedom
    # only fit submodels while the degrees of freedom does not become 
    # larger than the requested maximum
    if(df[sMax + 1] > dfMax) {
      keep <- which(df <= dfMax)
      s <- s[keep]
      sMax <- s[length(s)]
      firstActive <- firstActive[keep]
      df <- df[keep]
    }
    ## adjust for requested minimum step
    if(sMin > sMax) sMin <- sMax
    if(sMin > 0) {
      keep <- which(s >= sMin)
      s <- s[keep]
      df <- df[keep]
    }
    # groupwise sequenced variables (including intercept)
    sequenced <- c(1, unlist(assign[firstActive], use.names=FALSE) + 1)
  } else {
    # the default is to fit models as long as there are twice as many 
    # observations as predictors
    if(is.na(sMax)) sMax <- floor(n/2)
    if(sMax > length(active)) sMax <- length(active)
    if(sMin > sMax) sMin <- sMax
    s <- sMin:sMax
    # compute degrees of freedom of the submodels along sequence
    df <- s + 1  # account for intercept
    # sequenced variables (including intercept)
    sequenced <- c(1, active[seq_len(sMax)] + 1)
  }
  if(length(s) == 1) crit <- "none"
  # define function to fit the submodels along the sequence
  if(robust) {
    if(useFormula) {
      fitFun <- function(k) {
        x <- x[, sequenced[seq_len(k)], drop=FALSE]
        callRegFun(y ~ x - 1, fun=regFun, args=regArgs)
      }
    } else {
      fitFun <- function(k) {
        x <- x[, sequenced[seq_len(k)], drop=FALSE]
        callRegFun(x, y, fun=regFun, args=regArgs)
      }
    }
  } else {
    fitFun <- function(k) lm.fit(x[, sequenced[seq_len(k)], drop=FALSE], y)
  }
  # fit submodels
  # number of variables to use is one less than degrees of freedom
  if(is.null(cl)) models <- lapply(df, fitFun) 
  else models <- parLapply(cl, df, fitFun)
  # construct matrix of coefficents
  coef <- matrix(0, nrow=ncol(x), ncol=length(s), 
                 dimnames=list(colnames(x), s))
  for(k in seq_len(ncol(coef))) {
    coef[sequenced[seq_len(df[k])], k] <- coef(models[[k]])
  }
  # extract fitted values and residuals along sequence
  fitted <- sapply(models, fitted)
  residuals <- sapply(models, residuals)
  colnames(fitted) <- colnames(residuals) <- s
  # construct return object
  out <- list(active=active, s=s, coefficients=coef, fitted.values=fitted, 
              residuals=residuals, df=df, robust=robust)
  class(out) <- "seqModel"
  # add robust scale estimates if applicable
  if(robust) out$scale <- sapply(models, getScale)
  # compute optimality criterion along sequence
  if(crit == "BIC") out$crit <- bicSelect(out)
  # return results
  out
}
