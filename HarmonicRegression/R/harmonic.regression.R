calculate.amp.phi <- function(a.cos, b.sin) {
  amp <- unname(sqrt(a.cos^2 + b.sin^2))
  phi <- unname(atan2(b.sin, a.cos) %% (2*pi))
  return(cbind(amp=amp, phi=phi))
}

calculate.ci.amp.phi <- function(amp, a.cos, b.sin, fit.res.ssr, ssx, df) {
  if (det(ssx) == 0 | kappa(ssx) > 0.5/.Machine$double.eps)
    return(cbind(amp=NA, phi=NA))
  ssxinv <- solve(ssx)
  ssxinvab <- ssxinv[2:3, 2:3]
  Ja <- cbind(ifelse(amp > 0, a.cos/amp, 0), 
              ifelse(amp > 0, b.sin/amp, 0))
  Jp <- cbind(ifelse(amp > 0, -b.sin/amp^2, 0), 
              ifelse(amp > 0, a.cos/amp^2, 0))
  var.a <- apply(Ja, 1, function(x) x %*% ssxinvab %*% x) * fit.res.ssr/df
  var.p <- apply(Jp, 1, function(x) x %*% ssxinvab %*% x) * fit.res.ssr/df
  ##ciquant <- qnorm(0.025, lower.tail=FALSE)
  ciquant <- 2
  ci.amp <- sqrt(var.a)*ciquant
  ci.phi <- sqrt(var.p)*ciquant
  return(cbind(amp=ci.amp, phi=ifelse(ci.phi < pi, ci.phi, pi)))
}

harmonic.f.test <- function(rest.ssr, unrest.ssr, df) {
  ## f-statistic and pvalues
  fstat <- ((rest.ssr - unrest.ssr)/2) / (unrest.ssr/df)
  pval <- pf(fstat, 2, df, lower.tail=FALSE)
  return(pval)
}
    
harmonic.covar.matrix <- function(inputtime, Tau, mean.center=FALSE) {
  ## covariance matrix of independent variables
  X <- cbind(rep(1, length(inputtime)),
             cos(2*pi/Tau*inputtime), sin(2*pi/Tau*inputtime))
  if (mean.center)
    X <- sweep(X, 2, colMeans(X))
  ssxfull <- zapsmall(t(X) %*% X)
  return(ssxfull)
}
  
# normalize time series matrix (no NAs) -----------------------------------
normalize.ts.matrix <- function(inputts, inputtime, 
                                norm.pol, norm.pol.degree) {
  if (norm.pol) {
    if (length(inputtime) < (1 + norm.pol.degree))
      stop(paste("Normalization polynomial degree too high for the number",
                 "of time points"))
    trendfit <- lm(inputts ~ poly(inputtime, norm.pol.degree, raw=T))
    trend.ts <- fitted(trendfit)
    trend.coef <- coef(trendfit)
    if(any(zapsmall(trend.ts) == 0))
      stop(paste("Normalization using polynomial failed (zero-crossing);", 
                 "try other normalization settings"))
    return(list(norm.ts = inputts/trend.ts, norm.w = trend.coef, 
                norm.vals = trend.ts))
  }
  else {
    tsmeans <- colMeans(inputts)
    if(any(zapsmall(tsmeans) == 0))
      stop(paste("Some time series have zero means, normalization failed.",
                 "Try other normalization settings, or leave these data out."))
    return(list(norm.ts = scale(inputts, FALSE, tsmeans), norm.w = tsmeans))
  }
}


# harmonic regression matrix (no NAs) -------------------------------------
harmonic.regression.matrix <- function(inputts, inputtime, Tau) {
  
  ## check time series length
  if (length(inputtime) < 3)
    stop(paste("These time series are too short for a meaningful analysis.",
               "At least 3 time points are needed."))
  
  ## matrix fit of the unrestricted model (harmonic regression)
  inputts.fit <- lm(inputts ~ 1 + cos(2*pi/Tau*inputtime) + 
                      sin(2*pi/Tau*inputtime), x=TRUE)
  
  ## refrain from parameter estimation if the design matrix is bad
  sing.det <- zapsmall(inputts.fit$x)
  if (det(t(sing.det) %*% sing.det) == 0)
    stop(paste("The time points are so unfortunately spaced that a phase and",
               "amplitude determination is impossible."))
  
  ## fitted values, possibly coerce to matrix
  fit.vals <- as.matrix(fitted(inputts.fit))
  
  ## coefficients, amplitudes, phases
  coeffs <- t(coef(inputts.fit))
  pars <- as.data.frame(calculate.amp.phi(coeffs[, 2], coeffs[, 3]))
  if (!any(duplicated(colnames(inputts))))
    rownames(pars) <- colnames(inputts)
  ## covariance matrix of independent variables
  ssx <- harmonic.covar.matrix(inputtime, Tau)
  
  ## if more than 3 time points are available, confidence intervals and p-values
  ## can be computed.
  if (length(inputtime) == 3) {
    
    pvals <- NA
    ci <- NA
    fit.res.ssr <- NA
    
  } else {
  
    ## sum squared residual of the restricted and unrestricted models
    inputts.ssr <- apply(inputts, 2, var)*(length(inputtime) - 1)
    
    if (is.matrix(residuals(inputts.fit)))
      fit.res.ssr <- colSums(residuals(inputts.fit)^2)
    ## handle also the case with one single input ts
    else
      fit.res.ssr <- sum(residuals(inputts.fit)^2)
      
    
    ## f-statistic and pvalues
    pvals <- harmonic.f.test(inputts.ssr, fit.res.ssr, length(inputtime) - 3)
    
    ## distance to upper confidence interval limit
    if (abs(det(ssx)) > 0) {
      ci <- as.data.frame(calculate.ci.amp.phi(pars$amp, 
                                               coeffs[, 2], coeffs[, 3],
                                               fit.res.ssr, ssx, 
                                               length(inputtime) - 3))
      if (!any(duplicated(colnames(inputts))))
        rownames(ci) <- colnames(inputts)
    } else
      ci <- NA
    
  }
  
  ## return values
  return(list(fit.vals=fit.vals,
              pars=pars, pvals=pvals, qvals=p.adjust(pvals, method="BH"), 
              ci=ci, coeffs=coeffs[, 2:3], ssr=fit.res.ssr, 
  						df=(length(inputtime) - 3), ssx=ssx))
  
}


harmonic.regression.matrix.trend <- function(inputts, inputtime, Tau,
                                             trend.degree) {
  
  ## check for enough degrees of freedom
  df <- length(inputtime) - (trend.degree + 3)
  if (df < 0)
    stop(paste("Too few time points.  Unable to continue.  Try a lower setting",
               "for trend.degree."))
  
  ## matrix fit of the restricted model (polynomial)
  rest.fit <- lm(inputts ~ poly(inputtime, trend.degree, raw=T))
  
  ## matrix fit of the unrestricted model (harmonic regression)
  unrest.fit <- lm(inputts ~ cos(2*pi/Tau*inputtime) + 
                     sin(2*pi/Tau*inputtime) + poly(inputtime, trend.degree,
                                                    raw=T),
                   x=TRUE)
  ## possibly coerce fitted values to matrix
  fit.vals <- as.matrix(fitted(unrest.fit))
  
  ## refrain from parameter estimation if the design matrix is bad
  sing.det <- zapsmall(unrest.fit$x)
  if (det(t(sing.det) %*% sing.det) == 0)
    stop(paste("The time points are so unfortunately spaced that a phase and",
               "amplitude determination is impossible."))

  ## coefficients, amplitudes, phases
  coeffs <- t(coef(unrest.fit))
  pars <- as.data.frame(calculate.amp.phi(coeffs[, 2], coeffs[, 3]))
  if (!any(duplicated(colnames(inputts))))
    rownames(pars) <- colnames(inputts)
  ## covariance matrix of independent variables
  ssx <- harmonic.covar.matrix(inputtime, Tau)
  
  if (df == 0) {
    
    pvals <- NA
    ci <- NA
    unrest.ssr <- NA
    
  } else {
    
    ## sum squared residual of the restricted and unrestricted models
    if (is.matrix(residuals(rest.fit))) {
      rest.ssr <- colSums(residuals(rest.fit)^2)
      unrest.ssr <- colSums(residuals(unrest.fit)^2)
    } else {
      rest.ssr <- sum(residuals(rest.fit)^2)
      unrest.ssr <- sum(residuals(unrest.fit)^2)
    }
    
    ## f-statistic and pvalues
    pvals <- harmonic.f.test(rest.ssr, unrest.ssr, df)
    
    ## distance to upper confidence interval limit 
    if (abs(det(ssx)) > 0) {
      ci <- as.data.frame(calculate.ci.amp.phi(pars$amp, 
                                               coeffs[, 2], coeffs[, 3],
                                               unrest.ssr, ssx, df))
      if (!any(duplicated(colnames(inputts))))
        rownames(ci) <- colnames(inputts)
    } else
      ci <- NA
    
  }
  
  
  ## return values
  return(list(fit.vals=fit.vals,
              pars=pars, pvals=pvals, qvals=p.adjust(pvals, method="BH"), 
              ci=ci, coeffs=coeffs, ssr=unrest.ssr, df=df, ssx=ssx))
  
}


# normalize time series vector (with NAs) ---------------------------------
normalize.one.ts <- function(inputts, inputtime, 
                             norm.pol, norm.pol.degree) {
  n.non.na <- length(which(!is.na(inputts)))
  if (norm.pol) {
    if (n.non.na < (1 + norm.pol.degree))
      return(list(norm.ts = rep(NA, length(inputtime)), 
                  norm.w  = rep(NA, 1 + norm.pol.degree),
                  norm.vals = rep(NA, length(inputtime))))
    trendfit <- lm(inputts ~ poly(inputtime, norm.pol.degree, raw=T),
                   na.action=na.exclude)
    trend.ts <- fitted(trendfit)
    trend.coef <- coef(trendfit)
    if(any(zapsmall(trend.ts) == 0, na.rm=T)) {
      warning(paste("Zero-crossing of normalization polynomial.", 
                    "One of the time series is left out the fitting procedure"))
      return(list(norm.ts = rep(NA, length(inputts)), norm.w = trend.coef,
                  norm.vals = trend.ts))
    }
    return(list(norm.ts = inputts/trend.ts, norm.w = trend.coef,
                norm.vals = trend.ts))
  }
  else {
    if (n.non.na == 0)
      return(list(norm.ts=inputts, norm.w=NA))
    tsmean <- mean(inputts, na.rm=TRUE)
    if(tsmean < 10^-getOption("digits")) {
      warning(paste("Mean value zero of a time series.  It is left out of the",
                    "fitting procedure"))
      return(list(norm.ts = rep(NA, length(inputts)), norm.w = trend.ts))
    }
    return(list(norm.ts = inputts/tsmean, norm.w = tsmean))
  }
}


# harmonic regression one time series (with NAs) --------------------------
fit.one.harmonic <- function(inputts, inputtime, Tau) {
  
  n.non.na <- length(which(!is.na(inputts)))
  if (n.non.na < 3)
    return(list(pars = c(amp=NA, phi=NA),
                coeffs = rep(NA, 3), ci = c(amp=NA, phi=NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, df = NA, pval = NA))
  
  ## sum squared residual of the restricted model (mean centering)
  inputts.ssr <- var(inputts, na.rm=TRUE)*(n.non.na - 1)
  
  ## harmonic regression
  inputts.fit <- lm(inputts ~ (1 + cos(2*pi/Tau*inputtime) + 
                                 sin(2*pi/Tau*inputtime)),
                    na.action=na.exclude, x=TRUE) 
  
  ## refrain from parameter estimation if the design matrix is bad
  sing.det <- zapsmall(inputts.fit$x)
  if (det(t(sing.det) %*% sing.det) == 0)
    return(list(pars = c(amp=NA, phi=NA),
                coeffs = rep(NA, 3), ci = c(amp=NA, phi=NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, df = NA, pval = NA))
  
  
  fit.vals = fitted(inputts.fit)

  t.non.na <- inputtime[!is.na(inputts)] 
  ssx <- harmonic.covar.matrix(t.non.na, Tau)
 
  coeffs <- coef(inputts.fit)
  pars <- calculate.amp.phi(coeffs[2], coeffs[3])
  
  if (n.non.na == 3) {

    pval <- NA
    ci <- c(amp=NA, phi=NA)
    fit.res.ssr <- NA
  
  } else {
    
    ## sum squared residual of the unrestricted model; fitted values
    fit.res.ssr <- sum(residuals(inputts.fit)^2, na.rm=TRUE)
    
    ## f-statistic and pvalues
    pval <- harmonic.f.test(inputts.ssr, fit.res.ssr, n.non.na - 3)
    
    ## distance to upper confidence interval limit 
    ci <- calculate.ci.amp.phi(pars[, "amp"], coeffs[2], coeffs[3],
                               fit.res.ssr, ssx, n.non.na - 3)
    
  }
  
  return(list(pars = pars,
              coeffs = coeffs, ci = ci,
              fit.vals = fit.vals,
              ssr  = fit.res.ssr, df=(n.non.na - 3), ssx = ssx,
              pval = pval))

}


fit.one.harmonic.trend <- function(inputts, inputtime, Tau,
                                   trend.degree) {
  
  n.non.na <- length(which(!is.na(inputts)))
  ## check for enough degrees of freedom
  df <- n.non.na - (trend.degree + 3)
  if (df < 0)
    return(list(pars = c(amp=NA, phi=NA),
                coeffs = rep(NA, trend.degree + 3), ci = c(amp=NA, phi=NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, df = NA, pval = NA))
  
  ## fit of the restricted model (polynomial)
  rest.fit <- lm(inputts ~ poly(inputtime, trend.degree, raw=T))

  ## matrix fit of the unrestricted model (harmonic regression)
  unrest.fit <- lm(inputts ~ cos(2*pi/Tau*inputtime) + 
                     sin(2*pi/Tau*inputtime) + poly(inputtime, trend.degree,
                                                    raw=T),
                   na.action=na.exclude, x=TRUE)

  ## refrain from parameter estimation if the design matrix is bad
  sing.det <- zapsmall(unrest.fit$x)
  if (det(t(sing.det) %*% sing.det) == 0)
    return(list(pars = c(amp=NA, phi=NA),
                coeffs = rep(NA, trend.degree + 3), ci = c(amp=NA, phi=NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, df = NA, pval = NA))
  
  fit.vals <- fitted(unrest.fit)
  
  t.non.na <- inputtime[!is.na(inputts)] 
  ssx <- harmonic.covar.matrix(t.non.na, Tau)
  
  ## coefficients, amplitudes, phases
  coeffs <- coef(unrest.fit)
  pars <- calculate.amp.phi(coeffs[2], coeffs[3])
  
  if (df == 0) {

    pval <- NA
    ci <- c(amp=NA, phi=NA)
    unrest.ssr <- NA
  
  } else {
    
    ## sum squared residual of the restricted and unrestricted models
    rest.ssr <- sum(residuals(rest.fit)^2, na.rm=TRUE)
    unrest.ssr <- sum(residuals(unrest.fit)^2, na.rm=TRUE)

    ## f-statistic and pvalues
    pval <- harmonic.f.test(rest.ssr, unrest.ssr, df)
    
    ## distance to upper confidence interval limit (Halberg 1967)
    ci <- calculate.ci.amp.phi(pars[, "amp"], coeffs[2], coeffs[3],
                               unrest.ssr, ssx, df)
    
  }    
  
  return(list(pars = pars,
              coeffs = coef(unrest.fit), ci = ci,
              fit.vals = fit.vals,
              ssr  = unrest.ssr, df = df, ssx = ssx,
              pval = pval))
}


# harmonic regression one time series at a time, with NAs -----------------
harmonic.regression.nas <- function(inputts, inputtime, Tau) {
  
  if (length(inputtime) < 3)
    stop(paste("These time series are too short for a meaningful analysis.",
               "At least 3 time points are needed"))
  
  result.list <- apply(inputts, 2, fit.one.harmonic, inputtime, Tau)
  
  coeffs <- t(sapply(result.list, "[[", "coeffs"))
  fit.vals <- sapply(result.list, "[[", "fit.vals")
  ssr <- sapply(result.list, "[[", "ssr")
  df <- sapply(result.list, "[[", "df")
  ssx <- lapply(result.list, "[[", "ssx")
  pvals <- sapply(result.list, "[[", "pval")
  pars <- t(sapply(result.list, "[[", "pars"))
  colnames(pars) <- c("amp", "phi")
  if (!any(duplicated(colnames(inputts))))
    rownames(pars) <- colnames(inputts)
  pars <- as.data.frame(pars)
  ci <- t(sapply(result.list, "[[", "ci"))
  colnames(ci) <- c("amp", "phi")
  if (!any(duplicated(colnames(inputts))))
    rownames(ci) <- colnames(inputts)
  ci <- as.data.frame(ci)

  return(list(fit.vals=fit.vals,
              pars=pars, pvals=pvals, qvals=p.adjust(pvals, method="BH"),
              ci=ci, coeffs=coeffs, ssr=ssr, df=df, ssx=ssx))

}


harmonic.regression.nas.trend <- function(inputts, inputtime, Tau,
                                          trend.degree) {
  
  ## check for enough degrees of freedom
  if ((length(inputtime) - (trend.degree + 3)) < 0)
    stop(paste("Too few time points.  Unable to continue.  Try a lower", 
               "setting for trend.degree."))
  
  result.list <- apply(inputts, 2, fit.one.harmonic.trend, inputtime, Tau,
                       trend.degree)
  
  coeffs <- t(sapply(result.list, "[[", "coeffs"))
  fit.vals <- sapply(result.list, "[[", "fit.vals")
  ssr <- sapply(result.list, "[[", "ssr")
  df <- sapply(result.list, "[[", "df")
  ssx <- lapply(result.list, "[[", "ssx")
  pvals <- sapply(result.list, "[[", "pval")
  pars <- t(sapply(result.list, "[[", "pars"))
  colnames(pars) <- c("amp", "phi")
  if (!any(duplicated(colnames(inputts))))
    rownames(pars) <- colnames(inputts)
  pars <- as.data.frame(pars)
  ci <- t(sapply(result.list, "[[", "ci"))
  colnames(ci) <- c("amp", "phi")
  if (!any(duplicated(colnames(inputts))))
    rownames(ci) <- colnames(inputts)
  ci <- as.data.frame(ci)
  
  return(list(fit.vals=fit.vals,
              pars=pars, pvals=pvals, qvals=p.adjust(pvals, method="BH"),
              ci=ci, coeffs=coeffs, ssr=ssr, df=df, ssx=ssx))
  
}



# harmonic regression documentation ---------------------------------------
#' Harmonic Regression
#'
#' Estimates amplitudes and phases along with confidence intervals and p-values 
#' from a set of time series that may oscillate with a specified period. A 
#' model, per default \deqn{y = m + a cos(\omega t) + b sin(\omega t),} is 
#' fitted to the time series.  This model is equivalent to the model \deqn{m + c
#' cos(\omega t - \phi),} with amplitude \eqn{c = \sqrt(a^2 + b^2)} and phase 
#' \eqn{\phi = atan2(b, a)}. P-values for \eqn{c > 0} (more precisely: either 
#' \eqn{a} or \eqn{b > 0} ) are computed by an F-test.  Confidence intervals for
#' the amplitudes and phases are computed by a linear error propagation 
#' approximation.
#' 
#' The default setting is that the time series are normalized with their mean 
#' values.  Optionally a polynomial of degree 1 or more is first fitted to each
#' time series, whereupon the original time series are normalized by dividing
#' with the fitted values at each point, thus trends in a fold-change sense are
#' assumed.  Another option is trend elimination, in which case the same model
#' plus a polynomial: \eqn{y = m + a cos(\omega t) + b sin(\omega t) + et + ft^2
#' + ... } is fitted to the (possibly normalized) data.  In this case, returned
#' p-values still only concern the alternative \eqn{c > 0} as defined above.
#' 
#' Values returned include normalized time series (if normalization is 
#' performed), normalization weights (means or polynomial coefficients if
#' polynomial normalilzation is used), fitted normalized curves, fitted 
#' non-normalized curves, a data frame of amplitudes and phases (in radians), 
#' p-values according to an F-test (Halberg 1967), Benjamini-Hochberg adjusted 
#' p-values, a data frame of approximately 1.96 standard deviations for the 
#' amplitude and phase estimates, a matrix of coefficients a and b and possibly 
#' c,... , the sum square resuduals after the fit for each time series, and the 
#' covariance matrix for the three independent variables (\eqn{1}, 
#' \eqn{cos(\omega t)}, and \eqn{sin(\omega t)}).  The latter can be used in 
#' post-processing e.g. to obtain individual p-values for coefficients by 
#' t-tests.
#' 
#' @param inputts Matrix of time series.  Rows correspond to time points, 
#'   columns to samples.  If a vector is provided, it is coerced to a matrix.
#' @param inputtime Vector of the time points corresponding to the row in the 
#'   time series matrix.
#' @param Tau Scalar giving the oscillation period to estimate and test for.
#' @param normalize Boolean, set to \code{TRUE} if normalization is to be 
#'   performed (default).  Unless \code{norm.pol=TRUE}, normalization is
#'   performed by dividing with the mean.
#' @param norm.pol Boolean, set to \code{TRUE} if a polynomial should be fitted
#'   to each time series, and used for normalization.  In this case, each point
#'   in a time series will be divided by the value of the fitted polynomial.
#'   Defaults to \code{FALSE}.
#' @param norm.pol.degree Scalar indicating the polynomial degree for the
#'   normalization (ignored if \code{norm.pol=FALSE}).
#' @param trend.eliminate Boolean, set to \code{TRUE} if trend elimination is to
#'   be performed (see Details above).  Defaults to \code{FALSE}.
#' @param trend.degree Integer indicading the polynomial degree for the trend
#'   elimination, default is 1.  Ignored when \code{trend.eliminate=FALSE},
#'   which is the default.
#' @return A list containing: \tabular{ll}{
#'  \code{means} \tab Vector (if \code{norm.pol=FALSE}) or matrix (otherwise) of
#'  the means or coefficients of the fitted polynomial used for the
#'  normalization \cr
#'  \code{normts} \tab Matrix of mean-scaled or normalized-by-polynomial time 
#'  series, same dimensionality as \code{inputts} \cr
#'  \code{fit.vals} \tab Matrix of model fitted values to \code{inputts} \cr
#'  \code{norm.fit.vals} \tab Matrix of model fitted values to the normalized
#'  (trend eliminated or mean scaled) time series \cr
#'  \code{pars} \tab Data frame of estimated amplitudes and phases (in radians,
#'  between 0 and \eqn{2\pi}) \cr
#'  \code{pvals} \tab Vector of p-values according to an F-test of the model fit
#'  against a restricted model (mean-centering only) \cr
#'  \code{qvals} \tab Vector of Benjamini-Hochberg adjusted p-values \cr
#'  \code{ci} \tab Data frame of one-sided approximative 95\% (\eqn{2\sigma})
#'  confidence intervals for the estimated amplitudes and phases \cr
#'  \code{coeffs} \tab Matrix of estimated model parameters \eqn{a} and \eqn{b}
#'  \cr
#'  \code{ssr} \tab Vector of sum square residuals for the model fits \cr
#'  \code{df} \tab Scalar if \code{inputts} does not contain \code{NA}s and
#'  Vector otherwise, representing the degrees of freedom of the residual from
#'  the fit \cr
#'  \code{ssx} \tab Matrix (3 times 3, if \code{inputts} does not contain 
#'  \code{NA}s, a list of such matrices, one for each time series, otherwise) of
#'  covariances for the dependent variables corresponding to (\eqn{m}, \eqn{a
#'  cos(\omega t)}, and \eqn{b sin(\omega t)}, respecively) \cr
#' }
#' @export
#' @references Halberg F, Tong YL, Johnson EA: Circadian System Phase -- An 
#'   Aspect of Temporal Morphology; Procedures and Illustrative Examples. in:
#'   The Cellular Aspects of Biorhythms, Springer 1967.

# harmonic regression main function ---------------------------------------
harmonic.regression <- function(inputts, inputtime, Tau=24,
                                normalize=TRUE, norm.pol=FALSE, 
                                norm.pol.degree=1, trend.eliminate=FALSE,
                                trend.degree=1) {
  
  ## check that input data come as numerics
  if (!is.numeric(inputts) || !is.numeric(inputtime))
    stop("Input data (inputts, inputtime) must be numeric")
  
  ## try to coerce input to matrix, if needed
  if (is.vector(inputts))
    inputts <- as.matrix(inputts)
  
  ## check series lengths
  if (nrow(inputts) != length(inputtime))
    stop(paste("Length of time series (inputts):", nrow(inputts), " and time",
               "points (inputtime):", length(inputtime), "do not match."))
  
  if (norm.pol.degree < 1 || trend.degree < 1)
    stop(paste("Polynomials for normalization and trend elimination must both",
               "be of  degrees 1 or more. Aborting."))
  
  ## check if there are NAs
  if (any(is.na(inputts))) {

    ## if NAs; call the slow versions handling NAs separately for each time
    ## series
    if ("ts" %in% class(inputts))
      stop(paste("Time series object with NAs was supplied.  This is not yet",
                 "supported; please supply data containing NAs as a plain", 
                 "matrix.  Aborting."))
    
    if (normalize) {
      norm.ts.list <- apply(inputts, 2, normalize.one.ts,
                            inputtime, norm.pol, norm.pol.degree)
      norm.ts <- sapply(norm.ts.list, "[[", "norm.ts")
      norm.w <- sapply(norm.ts.list, "[[", "norm.w")
      if (trend.eliminate)
        results <- harmonic.regression.nas.trend(norm.ts, inputtime,
                                                 Tau, trend.degree)
      else
        results <- harmonic.regression.nas(norm.ts, inputtime, Tau)
      ## results$fit.vals are really fits to normalized time series
      ## fit.vals will be recalculated below
      results <- append(results, list(norm.fit.vals=results$fit.vals),
                        after=1)
      results <- append(results, list(means=norm.w),
                        after=0)
      results <- append(results, list(normts=norm.ts),
                        after=1)
      if (norm.pol) {
        norm.vals <- sapply(norm.ts.list, "[[", "norm.vals")
        results$fit.vals <- results$norm.fit.vals*norm.vals
      } else
        results$fit.vals <- sweep(results$norm.fit.vals, 2, norm.w, "*")
    } else {
      if (trend.eliminate)
        results <- harmonic.regression.nas.trend(inputts, inputtime, 
                                                 Tau, trend.degree)
      else
        results <- harmonic.regression.nas(inputts, inputtime, Tau)
      results <- append(results, list(means=colMeans(inputts, na.rm=TRUE)),
                        after=0)
    }
    
  } else {
    
    ## use fast vectorized versions
    
    if (normalize) {
      ## normalization
      norm.ts.list <- normalize.ts.matrix(inputts, inputtime, 
                                          norm.pol, norm.pol.degree)
      if (trend.eliminate)
        results <- harmonic.regression.matrix.trend(norm.ts.list$norm.ts,
                                                    inputtime, Tau, 
                                                    trend.degree)
      else
        results <- harmonic.regression.matrix(norm.ts.list$norm.ts, 
                                              inputtime, Tau)
      results <- append(results, list(norm.fit.vals=results$fit.vals),
                        after=1)
      results <- append(results, list(means=norm.ts.list$norm.w),
                        after=0)
      results <- append(results, list(normts=norm.ts.list$norm.ts),
                        after=1)
      ## compute non-normalized fitted values
      if (norm.pol) 
        results$fit.vals <- results$norm.fit.vals*norm.ts.list$norm.vals
      else
        results$fit.vals <- sweep(results$norm.fit.vals, 2, 
                                  norm.ts.list$norm.w, "*")
    } else {
      if (trend.eliminate)
        results <- harmonic.regression.matrix.trend(inputts,
                                                    inputtime, Tau, 
                                                    trend.degree)
      else
        results <- harmonic.regression.matrix(inputts, inputtime, Tau)
      results <- append(results, list(means=colMeans(inputts)),
                        after=0)
    }
  }
  
  return(results)

}



# Data sets documentation -------------------------------------------------

#' Menet et al. RNA-Seq Data
#' 
#' Quantification of circadian transcriptional activity in mouse liver.
#' 
#' @format a data frame with nascent RNA-seq data
#' @details The nascent-seq data are thought to reflect transcriptional
#'   activities. Raw data was collected from the supplementary material of the
#'   original publication cited below. In that study, samples were collected at
#'   6 different times of day in two biological replicates (fields ZT* in the
#'   data frame). The field 'mgi_symbol' refers to MGI gene names.
#' @source Jerome S Menet, Joseph Rodriguez, Katharine C Abruzzi, Michael 
#'   Rosbash. Nascent-Seq reveals novel features of mouse circadian 
#'   transcriptional regulation. eLife 1:e00011 (2012).
#' @examples 
#' data(rna.nasc)
#' nasc.t <- seq(0, 44, 4)
#' plot(nasc.t, rna.nasc["Arntl", -1], type="b")
#' @name rna.nasc
NULL


#' Menet et al. RNA-Seq Data
#' 
#' Quantification of circadian transcriptional activity in mouse liver.
#' 
#' @format a data frame with poly(A)+ RNA-seq data
#' @details The poly(A)+-seq data represent mature mRNA abundances.
#'   Raw data was collected from the supplementary material of the original 
#'   publication cited below. In that study, samples were collected at 6
#'   different times of day in two biological replicates (fields ZT* in the data
#'   frame). The field 'mgi_symbol' refers to MGI gene names.
#' @source Jerome S Menet, Joseph Rodriguez, Katharine C Abruzzi, Michael
#'   Rosbash. Nascent-Seq reveals novel features of mouse circadian
#'   transcriptional regulation. eLife 1:e00011 (2012).
#' @examples 
#' data(rna.polya)
#' polya.t <- seq(0, 44, 4)
#' plot(polya.t, rna.polya["Arntl", -1], type="b")
#' @name rna.polya
NULL









