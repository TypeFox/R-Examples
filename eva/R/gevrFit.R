#' Parameter estimation for the GEVr distribution model
#'
#' This function provides maximum likelihood estimation for the GEVr model, with the option of probability weighted moment and maximum product
#' spacing estimation for block maxima (GEV1) data. It also allows generalized linear modeling of the parameters.
#' @param data Data should be a matrix from the GEVr distribution.
#' @param method Method of estimation - maximum likelihood (mle), maximum product spacings (mps), and probability weighted moments (pwm). Uses mle by default.
#' For \eqn{r > 1}, only mle can be used.
#' @param information Whether standard errors should be calculated via observed or expected (default) information. For probability weighted moments,
#' only expected information will be used if possible. In the case with covariates, only observed information is available.
#' @param locvars,scalevars,shapevars A dataframe of covariates to use for modeling of the each parameter. Parameter
#' intercepts are automatically handled by the function. Defaults to NULL for the stationary model.
#' @param locform,scaleform,shapeform An object of class `formula' (or one that can be coerced into that class), specifying the model
#' of each parameter. By default, assumes stationary (intercept only) model. See details.
#' @param loclink,scalelink,shapelink A link function specifying the relationship between the covariates and each parameter. Defaults to the identity function. For
#' the stationary model, only the identity link should be used.
#' @param gumbel Whether to fit the Gumbel (type I) extreme value distribution (i.e. shape parameter equals zero). Defaults to FALSE.
#' @param start Option to provide a set of starting parameters to optim; a vector of location, scale, and shape, in that order. Otherwise, the routine attempts
#' to find good starting parameters. See details.
#' @param opt Optimization method to use with optim.
#' @param maxit Number of iterations to use in optimization, passed to optim. Defaults to 10,000.
#' @param ... Additional arguments to pass to optim.
#' @examples
#' set.seed(7)
#' x1 <- rgevr(500, 1, loc = 0.5, scale = 1, shape = 0.3)
#' result1 <- gevrFit(x1, method = "mps")
#'
#' ## A linear trend in the location and scale parameter
#' n <- 100
#' r <- 10
#' x2 <- rgevr(n, r, loc = 100 + 1:n / 50,  scale = 1 + 1:n / 300, shape = 0)
#'
#' covs <- as.data.frame(seq(1, n, 1))
#' names(covs) <- c("Trend1")
#' ## Create some unrelated covariates
#' covs$Trend2 <- rnorm(n)
#' covs$Trend3 <- 30 * runif(n)
#' result2 <- gevrFit(data = x2, method = "mle", locvars = covs, locform = ~ Trend1 + Trend2*Trend3,
#' scalevars = covs, scaleform = ~ Trend1)
#'
#' ## Show summary of estimates
#' result2
#'
#' @return A list describing the fit, including parameter estimates and standard errors for the mle and mps methods. Returns as a class
#' object `gevrFit' to be used with diagnostic plots.
#' @details In the stationary case (no covariates), starting parameters for mle and mps estimation are the probability weighted moment estimates.
#' In the case where covariates are used, the starting intercept parameters are the probability weighted moment estimates from the stationary case
#' and the parameters based on covariates are initially set to zero. For non-stationary parameters, the first reported estimate refers to the
#' intercept term. Covariates are centered and scaled automatically to speed up optimization, and then transformed back to original scale. \cr
#' Formulas for generalized linear modeling of the parameters should be given in the form `~ var1 + var2 + \eqn{\cdots}'. Essentially, specification
#' here is the same as would be if using function `lm' for only the right hand side of the equation. Interactions, polynomials, etc. can be
#' handled as in the `formula' class. \cr
#' Intercept terms are automatically handled by the function. By default, the link functions are the identity function and the covariate dependent
#' scale parameter estimates are forced to be positive. For some link function \eqn{f(\cdot)} and for example, scale
#' parameter \eqn{\sigma}, the link is written as \eqn{\sigma = f(\sigma_1 x_1 + \sigma_2 x_2 + \cdots + \sigma_k x_k)}. \cr
#' Maximum likelihood estimation can be used in all cases. Probability weighted moment estimation can only be used if \eqn{r = 1} and data is
#' assumed to be stationary. Maximum product spacings estimation can be used in the non-stationary case, but only if \eqn{r = 1}.
#'
#' @importFrom Matrix rankMatrix
#' @export
gevrFit <- function(data, method = c("mle", "mps", "pwm"), information = c("expected", "observed"), locvars = NULL, scalevars = NULL,
                     shapevars = NULL, locform = ~ 1, scaleform = ~ 1, shapeform = ~ 1, loclink = identity, scalelink = identity,
                     shapelink = identity, gumbel = FALSE, start = NULL, opt = "Nelder-Mead", maxit = 10000, ...) {
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  method <- match.arg(method)
  information <- match.arg(information)
  if((gumbel) & ((!is.null(shapevars)) | (shapeform != ~ 1)))
    stop("Cannot specify covariates in shape parameter for Gumbel distribution!")
  if(!is.null(locvars))
    if(nrow(locvars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if(!is.null(scalevars))
    if(nrow(scalevars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if(!is.null(shapevars))
    if(nrow(shapevars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if(((locform != ~ 1) & is.null(locvars)) | ((scaleform != ~ 1) & is.null(scalevars)) | ((shapeform != ~ 1) & is.null(shapevars)))
    stop("Need to specify covariates!")
  if(R > 1 & (method == "mps" | method == "pwm"))
    stop("If R > 1, MLE must be used")
  if((!is.null(locvars) | !is.null(scalevars) | !is.null(shapevars)) & method == "pwm")
    stop("Probability weighted moments can only be fitted for stationary data")

  if(locform == ~ 1)
    locvars <- as.data.frame(rep(1, n))

  if(scaleform == ~ 1)
    scalevars <- as.data.frame(rep(1, n))

  if(shapeform == ~ 1)
    shapevars <- as.data.frame(rep(1, n))

  locvars.model <- model.matrix(locform, data = locvars)
  locnames <- colnames(locvars.model)
  loccheck <- adjScale(locvars.model)
  if((rankMatrix(locvars.model)[1] < ncol(locvars.model)) | (rankMatrix(locvars.model)[1] > nrow(locvars.model)))
    stop("Location design matrix is singular")
  locvars.model <- loccheck$mat
  loctrans1 <- loccheck$adjmeans
  loctrans2 <- loccheck$adjvars

  scalevars.model <- model.matrix(scaleform, data = scalevars)
  scalenames <- colnames(scalevars.model)
  scalecheck <- adjScale(scalevars.model)
  if((rankMatrix(scalevars.model)[1] < ncol(scalevars.model)) | (rankMatrix(scalevars.model)[1] > nrow(scalevars.model)))
    stop("Scale design matrix is singular")
  scalevars.model <- scalecheck$mat
  scaletrans1 <- scalecheck$adjmeans
  scaletrans2 <- scalecheck$adjvars

  shapevars.model <- model.matrix(shapeform, data = shapevars)
  shapenames <- colnames(shapevars.model)
  shapecheck <- adjScale(shapevars.model)
  if((rankMatrix(shapevars.model)[1] < ncol(shapevars.model)) | (rankMatrix(shapevars.model)[1] > nrow(shapevars.model)))
    stop("Shape design matrix is singular")
  shapevars.model <- shapecheck$mat
  shapetrans1 <- shapecheck$adjmeans
  shapetrans2 <- shapecheck$adjvars

  trans1 <- c(loctrans1, scaletrans1, shapetrans1)
  trans2 <- c(loctrans2, scaletrans2, shapetrans2)

  locvars.model.orig <- t((t(locvars.model) * loctrans2) + loctrans1)
  scalevars.model.orig <- t((t(scalevars.model) * scaletrans2) + scaletrans1)
  shapevars.model.orig <- t((t(shapevars.model) * shapetrans2) + shapetrans1)

  ## Probability Weighted Moments
  ## Also use this as the intial estimates for other methods
  x <- sort(as.vector(data[, 1]))
  solve_shape <- function(sh, mom0, mom1, mom2) {
    (3^sh - 1) / (2^sh - 1) - (3 * mom2 - mom0) / (2 * mom1 - mom0)
  }
  moments <- rep(0, 3)
  for(j in 1:n) {
    moments[1] <- moments[1] + (x[j] / n)
    moments[2] <- moments[2] + ((j - 1) / (n - 1)) * (x[j] / n)
    moments[3] <- moments[3] + (((j - 1) * (j - 2)) / ((n - 1) * (n - 2))) * (x[j] / n)
  }
  mom0 <- moments[1]
  mom1 <- moments[2]
  mom2 <- moments[3]

  if(!gumbel) {
    shape0 <- uniroot(solve_shape, interval = c(-5, +5), mom0 = mom0, mom1 = mom1, mom2 = mom2)$root
    scale0 <- ((2 * mom1 - mom0) * shape0) / (gamma(1 - shape0) * (2^shape0 - 1))
    loc0 <- mom0 + scale0 * (1 - gamma(1 - shape0)) / shape0
    theta0 <- c(loc0, scale0, shape0)
  } else {
    scale0 <- (mom0 - 2 * mom1) / - log(2)
    loc0 <- mom0 + scale0 * digamma(1)
    theta0 <- c(loc0, scale0)
  }

  if(is.null(start)) {
    locinit <- c(loc0, rep(0, ncol(locvars.model) - 1))
    scaleinit <- c(scale0, rep(0, ncol(scalevars.model) - 1))
    if(!gumbel) {
      shapeinit <- c(shape0, rep(0, ncol(shapevars.model) - 1))
      init <- c(locinit, scaleinit, shapeinit)
    } else {
      init <- c(locinit, scaleinit)
    }
  } else {
    init <- start
  }

  parnum <- c(ncol(locvars.model), ncol(scalevars.model), ncol(shapevars.model))
  if(gumbel) parnum[3] <- 0

  if(method == "pwm") {
    names(theta0) <- c("Location", "Scale", "Shape")[1:length(theta0)]
    out <- list(n = n, data = data, par.ests = theta0, par.ses = NA, varcov = NA,
                converged = NA, nllh.final = NA, R = R,
                stationary = TRUE, parnum = parnum,
                par.sum = theta0, gumbel = gumbel,
                covars = list(locvars.model.orig, scalevars.model.orig, shapevars.model.orig),
                links = list(loclink, scalelink, shapelink),
                method = method, information = information)
  }

  negloglik <- function(vars, locvars1, scalevars1, shapevars1) {
    loc <- vars[1:length(locinit)]
    scale <- vars[(length(locinit) + 1):(length(locinit) + length(scaleinit))]
    if(!gumbel) shape <- vars[(length(locinit) + length(scaleinit) + 1):length(vars)] else shape <- 0
    locmat <- t(loc * t(locvars1))
    scalemat <- t(scale * t(scalevars1))
    shapemat <- t(shape * t(shapevars1))
    locvec <- loclink(rowSums(locmat))
    scalevec <- scalelink(rowSums(scalemat))
    shapevec <- shapelink(rowSums(shapemat))
    w <- matrix(((data - locvec) / scalevec), ncol = R)
    z <- pmax(matrix(w * shapevec, ncol = R), -1)
    cond1 <- any(scalevec <= 0)
    cond2 <- min(1 + w * shapevec) <= 0
    log.density <- ifelse(shapevec == 0, rowSums(-log(pmax(scalevec, 0)) - w) - exp(-w[,R]),
                          rowSums(-log(pmax(scalevec, 0)) - ((1/shapevec) + 1) * log1p(z)) -
                          exp((-1/shapevec) * log1p(z[,R])))
    log.density[(is.nan(log.density) | is.infinite(log.density))] <- 0
    if(cond1 | cond2) {
      abs(sum(log.density)) + 1e6
    } else {
      - sum(log.density)
    }
  }

  mpsobj <- function(vars, locvars1, scalevars1, shapevars1) {
    loc <- vars[1:length(locinit)]
    scale <- vars[(length(locinit) + 1):(length(locinit) + length(scaleinit))]
    if(!gumbel) shape <- vars[(length(locinit) + length(scaleinit) + 1):length(vars)] else shape <- 0
    locmat <- t(loc * t(locvars1))
    scalemat <- t(scale * t(scalevars1))
    shapemat <- t(shape * t(shapevars1))
    locvec <- loclink(rowSums(locmat))
    scalevec <- scalelink(rowSums(scalemat))
    shapevec <- shapelink(rowSums(shapemat))
    w <- as.vector((data - locvec) / scalevec)
    z <- pmax(w * shapevec, -1)
    cond1 <- any(scalevec <= 0)
    cond2 <- min(1 + w * shapevec) <= 0
    cdf <- ifelse(shapevec == 0, exp(-exp(-w)), exp(-exp((-1/shapevec)*log1p(z))))
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
                 locvars1 = locvars.model, scalevars1 = scalevars.model, shapevars1 = shapevars.model)

    if(fit$convergence)
      warning("optimization may not have succeeded")

    loc.ests <- fit$par[1:length(locinit)] / loctrans2
    scale.ests <- fit$par[(length(locinit) + 1):(length(locinit) + length(scaleinit))] / scaletrans2
    if(!gumbel)
      shape.ests <- fit$par[(length(locinit) + length(scaleinit) + 1):length(fit$par)] / shapetrans2

    loc.ests <- ifelse(loccheck$truevars == 0, loc.ests - sum(loc.ests * loctrans1), loc.ests)
    scale.ests <- ifelse(scalecheck$truevars == 0, scale.ests - sum(scale.ests * scaletrans1), scale.ests)
    if(!gumbel)
      shape.ests <- ifelse(shapecheck$truevars == 0, shape.ests - sum(shape.ests * shapetrans1), shape.ests)

    if(!gumbel) par.ests <- c(loc.ests, scale.ests, shape.ests) else par.ests <- c(loc.ests, scale.ests)

    if((information == "observed") | (locform != ~ 1) | (scaleform != ~ 1) | (shapeform != ~ 1)) {
      varcov <- solve(optimHess(par.ests, objfun, locvars1 = locvars.model.orig,
                                scalevars1 = scalevars.model.orig, shapevars1 = shapevars.model.orig))
    } else {
      varcov <- gevrFisher(data, par.ests, gumbel)
    }
    par.ses <- sqrt(diag(varcov))

    if(!gumbel) {
      names(par.ests) <- c(paste('Location', colnames(locvars.model.orig), sep = ' '),
                           paste('Scale', colnames(scalevars.model.orig), sep = ' '),
                           paste('Shape', colnames(shapevars.model.orig), sep = ' '))
    } else {
      names(par.ests) <- c(paste('Location', colnames(locvars.model.orig), sep = ' '),
                           paste('Scale', colnames(scalevars.model.orig), sep = ' '))
    }

    names(par.ses) <- names(par.ests)

    par.sum <- data.frame(par.ests, par.ses, par.ests / par.ses, 2 * pnorm(abs(par.ests / par.ses), lower.tail = FALSE))
    colnames(par.sum) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    par.sum$codes <- ifelse(par.sum[, 4] < 0.001, '***',
                            ifelse(par.sum[, 4] < 0.01, '**',
                                   ifelse(par.sum[, 4] < 0.05, '*',
                                          ifelse(par.sum[, 4] < 0.1, '.', ' '))))
    colnames(par.sum) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "")

    if(method == "mle") {
      out <- list(n = n, data = data, par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                  converged = fit$convergence, nllh.final = fit$value, R = R,
                  stationary = ((locform == ~ 1) & (scaleform == ~ 1) & (shapeform == ~ 1)),
                  parnum = parnum, par.sum = par.sum, gumbel = gumbel,
                  covars = list(locvars.model.orig, scalevars.model.orig, shapevars.model.orig),
                  links = list(loclink, scalelink, shapelink),
                  method = method, information = information)
    } else {
      out <- list(n = n, data = as.matrix(data), par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                  converged = fit$convergence, moran = fit$value, R = R,
                  stationary = ((locform == ~ 1) & (scaleform == ~ 1) & (shapeform == ~ 1)),
                  parnum = parnum, par.sum = par.sum, gumbel = gumbel,
                  covars = list(locvars.model.orig, scalevars.model.orig, shapevars.model.orig),
                  links = list(loclink, scalelink, shapelink),
                  method = method, information = information)
    }

  }

  class(out) <- "gevrFit"
  out
}


## S3 functions for class gevrFit
#' @export
plot.gevrFit <- function(x, ...) {
  gevrDiag(x, ...)
}


#' @export
print.gevrFit <- function(x, ...) {
  cat("Summary of fit:\n")
  print(x$par.sum, digits = 5)
  if(x$method != "pwm")
    cat("---\nSignif. codes:  0 '***' 0.001 '*' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
}


#' @export
logLik.gevrFit <- function (object, ...) {
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

