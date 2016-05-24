findthresh <- function(data, ne) {
  data <- rev(sort(as.numeric(data)))
  data[ne] - min(min(abs(diff(data))[abs(diff(data)) > 0]), 1e-6)
}


#' Parameter estimation for the Generalized Pareto Distribution (GPD)
#'
#' Fits exceedances above a chosen threshold to the Generalized Pareto model. Various estimation procedures can be used,
#' including maximum likelihood, probability weighted moments, and maximum product spacing. It also allows
#' generalized linear modeling of the parameters.
#' @param data Data should be a numeric vector from the GPD.
#' @param threshold A threshold value or vector of the same length as the data.
#' @param nextremes Number of upper extremes to be used (either this or the threshold must be given, but not both).
#' @param npp Length of each period (typically year). Is used in return level estimation. Defaults to 365.
#' @param method Method of estimation - maximum likelihood (mle), maximum product spacing (mps), and
#' probability weighted moments (pwm). Uses mle by default. For pwm, only the stationary model can be fit.
#' @param information Whether standard errors should be calculated via observed or expected (default) information. For probability
#' weighted moments, only expected information will be used if possible. For non-stationary models, only observed
#' information is used.
#' @param scalevars,shapevars A dataframe of covariates to use for modeling of the each parameter. Parameter
#' intercepts are automatically handled by the function. Defaults to NULL for the stationary model.
#' @param scaleform,shapeform An object of class `formula' (or one that can be coerced into that class), specifying the model
#' of each parameter. By default, assumes stationary (intercept only) model. See details.
#' @param scalelink,shapelink A link function specifying the relationship between the covariates and each parameter. Defaults to
#' the identity function. For the stationary model, only the identity link should be used.
#' @param start Option to provide a set of starting parameters to optim; a vector of scale and shape, in that order. Otherwise,
#' the routine attempts to find good starting parameters. See details.
#' @param opt Optimization method to use with optim.
#' @param maxit Number of iterations to use in optimization, passed to optim. Defaults to 10,000.
#' @param ... Additional arguments to pass to optim.
#' @return A class object `gpdFit' describing the fit, including parameter estimates and standard errors.
#' @details The base code for finding probability weighted moments is taken from the R package evir. See citation.
#' In the stationary case (no covariates), starting parameters for mle and mps estimation are the probability weighted moment estimates.
#' In the case where covariates are used, the starting intercept parameters are the probability weighted moment estimates from the
#' stationary case and the parameters based on covariates are initially set to zero. For non-stationary parameters, the
#' first reported estimate refers to the intercept term. Covariates are centered and scaled automatically to speed up optimization,
#' and then transformed back to original scale. \cr
#' Formulas for generalized linear modeling of the parameters should be given in the form `~ var1 + var2 + \eqn{\cdots}'. Essentially,
#' specification here is the same as would be if using function `lm' for only the right hand side of the equation. Interactions,
#' polynomials, etc. can be handled as in the `formula' class. \cr
#' Intercept terms are automatically handled by the function. By default, the link functions are the identity function and
#' the covariate dependent scale parameter estimates are forced to be positive. For some link function \eqn{f(\cdot)} and for
#' example, scale parameter \eqn{\sigma}, the link is written as \eqn{\sigma = f(\sigma_1 x_1 + \sigma_2 x_2 + \cdots + \sigma_k x_k)}. \cr
#' Maximum likelihood estimation and maximum product spacing estimation can be used in all cases. Probability weighted moments
#' can only be used for stationary models.
#' @examples
#' ## Fit data using the three different estimation procedures
#' set.seed(7)
#' x <- rgpd(2000, loc = 0, scale = 2, shape = 0.2)
#' ## Set threshold at 4
#' mle_fit <- gpdFit(x, threshold = 4, method = "mle")
#' pwm_fit <- gpdFit(x, threshold = 4, method = "pwm")
#' mps_fit <- gpdFit(x, threshold = 4, method = "mps")
#' ## Look at the difference in parameter estimates and errors
#' mle_fit$par.ests
#' pwm_fit$par.ests
#' mps_fit$par.ests
#'
#' mle_fit$par.ses
#' pwm_fit$par.ses
#' mps_fit$par.ses
#'
#' ## A linear trend in the scale parameter
#' set.seed(7)
#' n <- 300
#' x2 <- rgpd(n, loc = 0, scale = 1 + 1:n / 200, shape = 0)
#'
#' covs <- as.data.frame(seq(1, n, 1))
#' names(covs) <- c("Trend1")
#'
#' result1 <- gpdFit(x2, threshold = 0, scalevars = covs, scaleform = ~ Trend1)
#'
#' ## Show summary of estimates
#' result1
#'
#' @references Pfaff, Bernhard, Alexander McNeil, and A. Stephenson. "evir: Extreme Values in R." R package version (2012): 1-7.
#' @importFrom Matrix rankMatrix
#' @export
gpdFit <- function(data, threshold = NA, nextremes = NA, npp = 365, method = c("mle", "mps", "pwm"),
                   information = c("expected", "observed"), scalevars = NULL, shapevars = NULL, scaleform = ~ 1,
                   shapeform = ~ 1, scalelink = identity, shapelink = identity, start = NULL, opt = "Nelder-Mead",
                   maxit = 10000, ...) {
  data <- as.numeric(data)
  n <- length(data)
  method <- match.arg(method)
  information <- match.arg(information)
  if(is.na(nextremes) && is.na(threshold))
    stop("Enter either a threshold or the number of upper extremes")
  if(!is.na(nextremes) && !is.na(threshold))
    stop("Enter EITHER a threshold or the number of upper extremes")
  if(length(threshold) > 1)
    if(length(threshold) != n)
      stop("Threshold vector must be the same length as the data")
  if(!is.null(scalevars))
    if(nrow(scalevars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if(!is.null(shapevars))
    if(nrow(shapevars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if(((scaleform != ~ 1) & is.null(scalevars)) | ((shapeform != ~ 1) & is.null(shapevars)))
    stop("Need to specify covariates!")
  if((!is.null(scalevars) | !is.null(shapevars)) & method == "pwm")
    stop("Probability weighted moments can only be fitted for stationary data")
  if(!is.na(nextremes))
    threshold <- findthresh(data, nextremes)

  exceedtrue <- data > threshold
  exceedances <- data[exceedtrue]
  excess <- (data - threshold)[exceedtrue]
  Nu <- length(excess)
  p.less.thresh <- 1 - Nu/n

  if(scaleform == ~ 1)
    scalevars <- as.data.frame(rep(1, Nu))

  if(shapeform == ~ 1)
    shapevars <- as.data.frame(rep(1, Nu))

  scalevars.model <- model.matrix(scaleform, data = scalevars[exceedtrue, , drop = FALSE])
  scalenames <- colnames(scalevars.model)
  scalecheck <- adjScale(scalevars.model)
  if((rankMatrix(scalevars.model)[1] < ncol(scalevars.model)) | (rankMatrix(scalevars.model)[1] > nrow(scalevars.model)))
    stop("Scale design matrix is singular")
  scalevars.model <- scalecheck$mat
  scaletrans1 <- scalecheck$adjmeans
  scaletrans2 <- scalecheck$adjvars

  shapevars.model <- model.matrix(shapeform, data = shapevars[exceedtrue, , drop = FALSE])
  shapenames <- colnames(shapevars.model)
  shapecheck <- adjScale(shapevars.model)
  if((rankMatrix(shapevars.model)[1] < ncol(shapevars.model)) | (rankMatrix(shapevars.model)[1] > nrow(shapevars.model)))
    stop("Shape design matrix is singular")
  shapevars.model <- shapecheck$mat
  shapetrans1 <- shapecheck$adjmeans
  shapetrans2 <- shapecheck$adjvars

  trans1 <- c(scaletrans1, shapetrans1)
  trans2 <- c(scaletrans2, shapetrans2)

  scalevars.model.orig <- t((t(scalevars.model) * scaletrans2) + scaletrans1)
  shapevars.model.orig <- t((t(shapevars.model) * shapetrans2) + shapetrans1)

  xbar <- mean(excess)
  a0 <- xbar
  gamma <- -0.35
  delta <- 0
  pvec <- ((1:Nu) + gamma)/(Nu + delta)
  a1 <- mean(sort(excess) * (1 - pvec))
  shape0 <- 2 - a0/(a0 - 2 * a1)
  scale0 <- (2 * a0 * a1)/(a0 - 2 * a1)

  if(is.null(start)) {
    scaleinit <- c(scale0, rep(0, ncol(scalevars.model) - 1))
    shapeinit <- c(shape0, rep(0, ncol(shapevars.model) - 1))
    init <- c(scaleinit, shapeinit)
  } else {
    init <- start
  }

  parnum <- c(ncol(scalevars.model), ncol(shapevars.model))

  if(method == "pwm") {
    denom <- Nu * (1 - 2 * shape0) * (3 - 2 * shape0)
    if(shape0 > 0.5) {
      denom <- NA
      warning("Asymptotic standard errors not available for",
              "PWM Method when shape > 0.5")
    }
    one <- (7 - 18 * shape0 + 11 * shape0^2 - 2 * shape0^3) * scale0^2
    two <- (1 - shape0) * (1 - shape0 + 2 * shape0^2) * (2 - shape0)^2
    cov <- scale0 * (2 - shape0) * (2 - 6 * shape0 + 7 * shape0^2 - 2 *
                                      shape0^3)
    varcov <- matrix(c(one, cov, cov, two), 2) / denom
    par.ses <- sqrt(diag(varcov))
    par.ests <- c(scale0, shape0)
  }

  negloglik <- function(vars, scalevars1, shapevars1) {
    scale <- vars[1:length(scaleinit)]
    shape <- vars[(length(scaleinit) + 1):length(vars)]
    scalemat <- t(scale * t(scalevars1))
    shapemat <- t(shape * t(shapevars1))
    scalevec <- scalelink(rowSums(scalemat))
    shapevec <- shapelink(rowSums(shapemat))
    w <- excess / scalevec
    cond1 <- any(scalevec <= 0)
    cond2 <- min(1 + w * shapevec) <= 0
    log.density <- -log(pmax(scalevec, 0)) -
      ifelse(shapevec == 0, w, ((1/shapevec) + 1) * log1p(pmax(w * shapevec, -1)))
    log.density[(is.nan(log.density) | is.infinite(log.density))] <- 0
    if(cond1 | cond2) {
      abs(sum(log.density)) + 1e6
    } else {
      - sum(log.density)
    }
  }

  mpsobj <- function(vars, scalevars1, shapevars1) {
    scale <- vars[1:length(scaleinit)]
    shape <- vars[(length(scaleinit) + 1):length(vars)]
    scalemat <- t(scale * t(scalevars1))
    shapemat <- t(shape * t(shapevars1))
    scalevec <- scalelink(rowSums(scalemat))
    shapevec <- shapelink(rowSums(shapemat))
    w <- excess / scalevec
    cond1 <- any(scalevec <= 0)
    cond2 <- min(1 + w * shapevec) <= 0
    cdf <- ifelse(shapevec == 0, 1 - exp(-w), 1 - exp((-1/shapevec)*log1p(pmax(w * shapevec, -1))))
    cdf[(is.nan(cdf) | is.infinite(cdf))] <- 0
    cdf <- sort(cdf)
    cdf <- c(0, cdf, 1)
    D <- diff(cdf)
    ## Check if any differences are zero due to rounding and adjust
    D <- ifelse(D == 0, .Machine$double.eps, D)
    if(cond1 | cond2) {
      abs(sum(log(D))) + 1e6
    } else {
      - sum(log(D))
    }
  }

  if(method == "mle")
    objfun <- negloglik
  if(method == "mps")
    objfun <- mpsobj

  if(method == "mle" | method == "mps") {
    fit <- optim(init, objfun, hessian = FALSE, method = opt, control = list(maxit = maxit, ...),
                 scalevars1 = scalevars.model, shapevars1 = shapevars.model)

    if(fit$convergence)
      warning("optimization may not have succeeded")

    scale.ests <- fit$par[1:length(scaleinit)] / scaletrans2
    shape.ests <- fit$par[(length(scaleinit) + 1):length(fit$par)] / shapetrans2

    scale.ests <- ifelse(scalecheck$truevars == 0, scale.ests - sum(scale.ests * scaletrans1), scale.ests)
    shape.ests <- ifelse(shapecheck$truevars == 0, shape.ests - sum(shape.ests * shapetrans1), shape.ests)

    par.ests <- c(scale.ests, shape.ests)

    if((information == "observed") | (scaleform != ~ 1) | (shapeform != ~ 1)) {
      varcov <- solve(optimHess(par.ests, objfun, scalevars1 = scalevars.model.orig,
                                shapevars1 = shapevars.model.orig))
    } else {
      varcov <- gpdFisher(Nu, par.ests)
    }
    par.ses <- sqrt(diag(varcov))
  }

  names(par.ests) <- c(paste('Scale', colnames(scalevars.model.orig), sep = ' '),
                       paste('Shape', colnames(shapevars.model.orig), sep = ' '))

  names(par.ses) <- names(par.ests)

  par.sum <- data.frame(par.ests, par.ses, par.ests / par.ses, 2 * pnorm(abs(par.ests / par.ses), lower.tail = FALSE))
  colnames(par.sum) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  par.sum$codes <- ifelse(par.sum[, 4] < 0.001, '***',
                          ifelse(par.sum[, 4] < 0.01, '**',
                                 ifelse(par.sum[, 4] < 0.05, '*',
                                        ifelse(par.sum[, 4] < 0.1, '.', ' '))))
  colnames(par.sum) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "")

  if(method == "mle") {
    out <- list(n = n, data = data, threshold = threshold, par.sum = par.sum,
                exceedances = exceedances, n.exceed = Nu, converged = fit$convergence,
                p.less.thresh = p.less.thresh, method = method, nllh.final = fit$value,
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                parnum = parnum, links = list(scalelink, shapelink),
                covars = list(scalevars.model.orig, shapevars.model.orig),
                stationary = ((scaleform == ~ 1) & (shapeform == ~ 1)),
                information = information, npp = npp, rate = 1 - p.less.thresh)
  }

  if(method == "mps") {
    out <- list(n = n, data = data, threshold = threshold, par.sum = par.sum,
                exceedances = exceedances, n.exceed = Nu, converged = fit$convergence,
                p.less.thresh = p.less.thresh, method = method, moran = fit$value,
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                parnum = parnum, links = list(scalelink, shapelink),
                covars = list(scalevars.model.orig, shapevars.model.orig),
                stationary = ((scaleform == ~ 1) & (shapeform == ~ 1)),
                information = information, npp = npp, rate = 1 - p.less.thresh)
  }

  if(method == "pwm") {
    out <- list(n = n, data = data, threshold = threshold, par.sum = par.sum,
                exceedances = exceedances, n.exceed = Nu,
                p.less.thresh = p.less.thresh, method = method,
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                parnum = parnum, links = list(scalelink, shapelink),
                covars = list(scalevars.model.orig, shapevars.model.orig),
                stationary = TRUE, information = "expected", npp = npp,
                rate = 1 - p.less.thresh)
  }

  class(out) <- "gpdFit"
  out
}


## S3 functions for class gpdFit
#' @export
plot.gpdFit <- function(x, ...) {
  gpdDiag(x, ...)
}


#' @export
print.gpdFit <- function(x, ...) {
  cat("Summary of fit:\n")
  print(x$par.sum, digits = 5)
  cat("---\nSignif. codes:  0 '***' 0.001 '*' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
}


#' @export
logLik.gpdFit <- function (object, ...) {
  if(!missing(...))
    warning("Extra arguments discarded")
  if(object$method != "mle")
    stop("Estimation method is not maximum likelihood")
  val <- - object$nllh.final
  attr(val, "nobs") <- object$n
  attr(val, "df") <- sum(object$parnum)
  class(val) <- "logLik"
  val
}
