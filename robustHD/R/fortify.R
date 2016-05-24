# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## prepare BIC values for plotting
# model .... object of class "bicSelect"
# data ..... data frame containing additional information such as step numbers 
#            or values of a tuning parameter
# select ... indicates columns of a BIC matrix to keep
fortify.bicSelect <- function(model, data = NULL, select = NULL, ...) {
  # extract BIC values and make sure they are in matrix form
  bic <- model$values
  if(!is.null(dim(bic)) && !is.null(select)) bic <- bic[, select, drop=FALSE]
  bic <- as.data.frame(bic)
  # define index of the submodels
  d <- dim(bic)
  # check data frames for BIC and additional information
  if(any(d == 0)) stop("BIC data has no rows or columns")
  if(is.null(data)) data <- bic[, 0, drop=FALSE]  # NULL data frame
  else if(!is.data.frame(data) || nrow(data) != d[1]) {
    stop(sprintf("'data' must be a data frame with %d rows", d[1]))
  }
  # combine all information into one data frame
  index <- seq_len(d[1])
  if(d[2] == 1) {
    bic <- cbind(index, data, bic)
    names(bic) <- c("index", names(data), "BIC")
  } else {
    # reshape BIC matrix and add additional information (works if data is NULL)
    bic <- mapply(function(val, nam) {
      cbind(fit=rep.int(nam, d[1]), index=index, data, BIC=val)
    }, val=bic, nam=names(bic), SIMPLIFY=FALSE, USE.NAMES=FALSE)
    bic <- do.call(rbind, bic)
    attr(bic, "facets") <- . ~ fit
  }
  # return data frame
  bic
}



#' Convert a sequence of regression models into a data frame for plotting
#' 
#' Supplement the fitted values and residuals of a sequence of regression 
#' models (such as robust least angle regression models or sparse least trimmed 
#' squares regression models) with other useful information for diagnostic 
#' plots.
#' 
#' @method fortify seqModel
#' @aliases fortify.rlars
#' 
#' @param model  the model fit to be converted.
#' @param data  currently ignored.
#' @param s  for the \code{"seqModel"} method, an integer vector giving 
#' the steps of the submodels to be converted (the default is to use the 
#' optimal submodel).  For the \code{"sparseLTS"} method, an integer vector 
#' giving the indices of the models to be converted (the default is to use the 
#' optimal model for each of the requested fits).
#' @param fit  a character string specifying which fit to convert.  Possible 
#' values are \code{"reweighted"} (the default) to convert the reweighted fit, 
#' \code{"raw"} to convert the raw fit, or \code{"both"} to convert both fits.
#' @param covArgs  a list of arguments to be passed to 
#' \code{\link[robustbase]{covMcd}} for computing robust Mahalanobis distances.
#' @param \dots  additional arguments to be passed to 
#' \code{\link[robustbase]{covMcd}} can be specified directly instead of via 
#' \code{covArgs}.
#' 
#' @return  A data frame containing the columns listed below, as well as 
#' additional information stored in the attributes \code{"qqLine"} (intercepts 
#' and slopes of the respective reference lines to be displayed in residual Q-Q 
#' plots), \code{"q"} (quantiles of the Mahalanobis distribution used as cutoff 
#' points for detecting leverage points), and \code{"facets"} (default faceting 
#' formula for the diagnostic plots).
#' @returnItem step  the steps (for the \code{"seqModel"} method) or indices 
#' (for the \code{"sparseLTS"} method) of the models (only returned if more 
#' than one model is requested).
#' @returnItem fit  the model fits (only returned if both the reweighted 
#' and raw fit are requested in the \code{"sparseLTS"} method).
#' @returnItem index  the indices of the observations.
#' @returnItem fitted  the fitted values.
#' @returnItem residual  the standardized residuals.
#' @returnItem theoretical  the corresponding theoretical quantiles from the 
#' standard normal distribution.
#' @returnItem qqd  the absolute distances from a reference line through the 
#' first and third sample and theoretical quartiles.
#' @returnItem rd  the robust Mahalanobis distances computed via the MCD (see 
#' \code{\link[robustbase]{covMcd}}).
#' @returnItem xyd  the pairwise maxima of the absolute values of the 
#' standardized residuals and the robust Mahalanobis distances, divided by the 
#' respective other outlier detection cutoff point.
#' @returnItem weight  the weights indicating regression outliers.
#' @returnItem leverage  logicals indicating leverage points (i.e., outliers in 
#' the predictor space).
#' @returnItem classification  a factor with levels \code{"outlier"} 
#' (regression outliers) and \code{"good"} (data points following the model).
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[ggplot2]{fortify}}, \code{\link{diagnosticPlot}}, 
#' \code{\link{sparseLTS}}, \code{\link{sparseLTS}}
#' 
#' @example inst/doc/examples/example-fortify.R
#' 
#' @keywords utilities
#' 
#' @import ggplot2
#' @export

fortify.seqModel <- function(model, data, s = NA, covArgs = list(...), ...) {
  ## initializations
  if(!model$robust) stop("not implemented yet")
  # check the scale estimate
  scale <- getComponent(model, "scale", s=s)
  if(any(scale <= 0)) stop("residual scale equal to 0")
  # check if model data is available to compute robust MCD distances
  terms <- delete.response(model$terms)  # extract terms for model matrix
  if(is.null(x <- model$x)) {
    x <- try(model.matrix(terms), silent=TRUE)
    if(inherits(x, "try-error")) {
      x <- NULL
      warning("model data not available")
    }
  }
  if(!is.null(x)) x <- removeIntercept(x)
  ## construct data frame with all information for plotting
  steps <- model$s
  if(length(steps) > 1) {
    ## check steps
    if(is.null(s)) s <- steps
    else if(isTRUE(is.na(s))) s <- getSOpt(model)  # defaults to optimal step
    else s <- checkSteps(s, sMin=steps[1], sMax=steps[length(steps)])
  } else s <- NA
  ## extract data for the requested steps
  if(length(s) > 1) {
    # extract the data from each requested step
    data <- lapply(s, fortifySeqModelStep, model=model, x=x, covArgs=covArgs)
    qql <- lapply(data, attr, which="qqLine")
    q <- lapply(data, attr, which="q")
    # combine data from the steps
    data <- cbind(step=rep.int(s, sapply(data, nrow)), do.call(rbind, data))
    qql <- cbind(step=rep.int(s, sapply(qql, nrow)), do.call(rbind, qql))
    q <- cbind(step=rep.int(s, sapply(q, nrow)), do.call(rbind, q))
    attr(data, "facets") <- ~step
    attr(data, "qqLine") <- qql
    attr(data, "q") <- q
  } else {
    # extract the data from the selected step
    data <- fortifySeqModelStep(s, model=model, x=x, covArgs=covArgs)
  }
  ## return data frame
  data
}


## workhorse function for a single step from a sequence of regression models
fortifySeqModelStep <- function(s, model, x = NULL, covArgs = list()) {
  ## extract fitted values
  fitted <- fitted(model, s=s)
  ## extract standardized residuals
  residuals <- residuals(model, s=s, standardized=TRUE)
  n <- length(residuals)  # number of observations
  ## compute 0/1 outlier weights
  wt <- as.integer(abs(residuals) <= qnorm(0.9875))
  ## compute theoretical quantiles and distances from Q-Q reference line
  theoretical <- qqNorm(residuals)
  qql <- qqLine(residuals)  # Q-Q reference line
  qqd <- abs(residuals - qql$intercept - qql$slope * theoretical)
  ## compute MCD distances using significant variables
  # extract predictor matrix
  ok <- (is.na(s) || s > 0) && !is.null(x)
  if(ok) {
    # extract coefficients
    coefficients <- removeIntercept(coef(model, s=s))
    significant <- which(coefficients != 0)
    p <- length(significant)
    if(p == 0) {
      ok <- FALSE
      warning("all coefficients equal to 0")
    }
  }
  if(ok) {
    # compute distances
    rd <- try({
      xs <- x[, significant, drop=FALSE]
      callCovFun <- getCallFun(covArgs)
      mcd <- callCovFun(xs, fun=covMcd, args=covArgs)
      sqrt(mahalanobis(xs, mcd$center, mcd$cov))
    }, silent=TRUE)
    if(inherits(rd, "try-error")) {
      ok <- FALSE
      warning("robust distances cannot be computed")
    }
  }
  if(!ok) rd <- rep.int(NA, n)
  # take maximum of the distances in the x- and y-space, divided by the 
  # respective other cutoff point
  q <- sqrt(qchisq(0.975, p))
  xyd <- pmax.int(abs(rd/2.5), abs(residuals/q))
  ## construct indicator variables for leverage points
  leverage <- rd > q
  ## classify data points
  class <- ifelse(wt == 0, "outlier", "good")
  class <- factor(class, levels=c("outlier", "good"))
  ## construct data frame
  data <- data.frame(index=seq_len(n), fitted=fitted, residual=residuals, 
                     theoretical=theoretical, qqd=qqd, rd=rd, xyd=xyd, 
                     weight=wt, leverage=leverage, classification=class)
  attr(data, "qqLine") <- as.data.frame(qql)
  attr(data, "q") <- data.frame(q=max(q, 2.5))
  ## return data frame
  data
}


#' @rdname fortify.seqModel
#' @method fortify sparseLTS
#' @export

fortify.sparseLTS <- function(model, data, s = NA, 
                              fit = c("reweighted", "raw", "both"), 
                              covArgs = list(...), ...) {
  ## initializations
  fit <- match.arg(fit)
  # check the scale estimate
  scale <- getComponent(model, "scale", s=s, fit=fit)
  if(any(scale <= 0)) stop("residual scale equal to 0")
  # check if model data is available to compute robust MCD distances
  terms <- delete.response(model$terms)  # extract terms for model matrix
  if(is.null(x <- model$x)) {
    x <- try(model.matrix(terms), silent=TRUE)
    if(inherits(x, "try-error")) {
      x <- NULL
      warning("model data not available")
    }
  }
  if(!is.null(x) && model$intercept) x <- removeIntercept(x)
  ## construct data frame with all information for plotting
  lambda <- model$lambda
  bothOpt <- FALSE
  if(length(lambda) > 1) {
    ## check steps
    steps <- seq_along(lambda)
    if(is.null(s)) s <- steps
    else if(isTRUE(is.na(s))) {
      s <- getSOpt(model, fit=fit)  # defaults to optimal step
      if(fit == "both") {
        s <- unique(s)
        bothOpt <- length(s) == 2
      }
    } else s <- checkSteps(s, sMin=1, sMax=length(steps))
  } else s <- NA
  ## extract data for the requested steps
  if(bothOpt) {
    # extract the data from the respecitve optimal lambda
    fits <- c("reweighted", "raw")
    ## recursive call for each fit
    reweighted <- fortifySparseLTSFit(model, s=s[1], fit="reweighted", x=x)
    raw <- fortifySparseLTSFit(model, s=s[2], fit="raw", x=x)
    ## combine data for Q-Q reference line
    qql <- data.frame(fit=factor(fits, levels=fits), 
                      rbind(attr(reweighted, "qqLine"), attr(raw, "qqLine")), 
                      row.names=NULL)
    ## combine data for cutoff chi-squared quantile
    q <- data.frame(fit=factor(fits, levels=fits), 
                    rbind(attr(reweighted, "q"), attr(raw, "q")), 
                    row.names=NULL)
    ## combine results
    n <- c(nrow(reweighted), nrow(raw))
    data <- data.frame(fit=rep.int(factor(fits, levels=fits), n), 
                       rbind(reweighted, raw), row.names=NULL)
    attr(data, "facets") <- . ~ fit
    attr(data, "qqLine") <- qql
    attr(data, "q") <- q
  } else if(length(s) > 1) {
    # extract the data from each requested step
    data <- lapply(s, fortifySparseLTSStep, model=model, fit=fit, 
                   x=x, covArgs=covArgs)
    qql <- lapply(data, attr, which="qqLine")
    q <- lapply(data, attr, which="q")
    # combine data from the steps
    data <- cbind(step=rep.int(s, sapply(data, nrow)), do.call(rbind, data))
    qql <- cbind(step=rep.int(s, sapply(qql, nrow)), do.call(rbind, qql))
    q <- cbind(step=rep.int(s, sapply(q, nrow)), do.call(rbind, q))
    attr(data, "facets") <- if(fit == "both") step ~ fit else ~step
    attr(data, "qqLine") <- qql
    attr(data, "q") <- q
  } else {
    # extract the data from the selected step
    data <- fortifySparseLTSStep(s, model=model, fit=fit, x=x, covArgs=covArgs)
  }
  ## return data frame
  data
}


## workhorse functions

# fortify a single sparse LTS step
fortifySparseLTSStep <- function(s, model, fit = "reweighted", 
                                 x = NULL, covArgs = list()) {
  ## construct data frame with all information for plotting
  if(fit == "both") {
    fits <- c("reweighted", "raw")
    ## call workhorse function for each fit
    reweighted <- fortifySparseLTSFit(model, s=s, fit="reweighted", 
                                      x=x, covArgs=covArgs)
    raw <- fortifySparseLTSFit(model, s=s, fit="raw", x=x, covArgs=covArgs)
    ## combine data for Q-Q reference line
    qql <- data.frame(fit=factor(fits, levels=fits), 
                      rbind(attr(reweighted, "qqLine"), attr(raw, "qqLine")), 
                      row.names=NULL)
    ## combine data for cutoff chi-squared quantile
    q <- data.frame(fit=factor(fits, levels=fits), 
                    rbind(attr(reweighted, "q"), attr(raw, "q")), 
                    row.names=NULL)
    ## combine results
    n <- c(nrow(reweighted), nrow(raw))
    data <- data.frame(fit=rep.int(factor(fits, levels=fits), n), 
                       rbind(reweighted, raw), row.names=NULL)
    attr(data, "facets") <- . ~ fit
    attr(data, "qqLine") <- qql
    attr(data, "q") <- q
  } else data <- fortifySparseLTSFit(model, s=s, fit=fit, x=x, covArgs=covArgs)
  ## return data
  data
}

# fortify a single sparse LTS fit
fortifySparseLTSFit <- function(model, s, fit = "reweighted", 
                                x = NULL, covArgs = list()) {
  ## extract fitted values
  fitted <- fitted(model, s=s, fit=fit)
  ## extract standardized residuals
  residuals <- residuals(model, s=s, fit=fit, standardized=TRUE)
  n <- length(residuals)  # number of observations
  ## extract outlier weights
  wt <- wt(model, s=s, fit=fit)
  ## compute theoretical quantiles and distances from Q-Q reference line
  theoretical <- qqNorm(residuals)
  qql <- qqLine(residuals)  # Q-Q reference line
  qqd <- abs(residuals - qql$intercept - qql$slope * theoretical)
  ## compute MCD distances using significant variables
  # extract predictor matrix
  ok <- !is.null(x)
  if(ok) {
    # extract coefficients
    coefficients <- coef(model, s=s, fit=fit)
    if(model$intercept) coefficients <- removeIntercept(coefficients)
    significant <- which(coefficients != 0)
    p <- length(significant)
    if(p == 0) {
      ok <- FALSE
      warning("all coefficients equal to 0")
    }
  }
  if(ok) {
    # adjust alpha since MCD computes subset size depending on n and p
    h <- model$quan
    n2 <- (n+p+1) %/% 2
    alpha <- pmin((h - 2*n2 + n) / (2 * (n - n2)), 1)
    # check fraction for subset size
    if(alpha < 0.5) {
      alpha <- 0.5
      warning(sprintf("cannot compute MCD with h = %d; using h = %d", 
                      model$quan, h.alpha.n(alpha, n, p)))
    }
    # compute distances
    rd <- try({
      xs <- x[, significant, drop=FALSE]
      covArgs$alpha <- alpha
      callCovFun <- getCallFun(covArgs)
      mcd <- callCovFun(xs, fun=covMcd, args=covArgs)
      if(fit == "reweighted") {
        center <- mcd$center
        cov <- mcd$cov
      } else {
        center <- mcd$raw.center
        cov <- mcd$raw.cov
      }
      sqrt(mahalanobis(xs, center, cov))
    }, silent=TRUE)
    if(inherits(rd, "try-error")) {
      ok <- FALSE
      warning("robust distances cannot be computed")
    }
  }
  if(!ok) rd <- rep.int(NA, n)
  # take maximum of the distances in the x- and y-space, divided by the 
  # respective other cutoff point
  q <- sqrt(qchisq(0.975, p))
  xyd <- pmax.int(abs(rd/2.5), abs(residuals/q))
  ## construct indicator variables for leverage points
  leverage <- rd > q
  ## classify data points
  class <- ifelse(wt == 0, "outlier", "good")
  class <- factor(class, levels=c("outlier", "good"))
  ## construct data frame
  data <- data.frame(index=seq_len(n), fitted=fitted, residual=residuals, 
                     theoretical=theoretical, qqd=qqd, rd=rd, xyd=xyd, 
                     weight=wt, leverage=leverage, classification=class)
  attr(data, "qqLine") <- as.data.frame(qql)
  attr(data, "q") <- data.frame(q=max(q, 2.5))
  ## return data frame
  data
}
