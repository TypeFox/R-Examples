############################
# S3 method for gmnl package
#############################

#' @rdname gmnl
#' @method print gmnl
#' @import stats
#' @export
print.gmnl <- function(x, digits = max(3, getOption("digits") - 3),
                          width = getOption("width"), ...){
  cat("\nCall:\n", deparse(x$call),"\n\n", sep = "")
  
  cat("\nCoefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}


#' @rdname gmnl
#' @method summary gmnl
#' @import stats
#' @export
summary.gmnl <- function(object,...){
  b <- object$coefficients
  std.err <- sqrt(diag(vcov(object)))
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  class(object) <- c("summary.gmnl", "gmnl")
  return(object)
}

#' @rdname gmnl
#' @method print summary.gmnl
#' @import stats
#' @export
print.summary.gmnl <- function(x, digits = max(3, getOption("digits") - 2),
                                  width = getOption("width"),
                                  ...){
  cat(paste("\nModel estimated on:", format(Sys.time(), "%a %b %d %X %Y"), "\n"))
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  
  cat("\nFrequencies of categories:\n")
  print(prop.table(x$freq), digits = digits)
  cat("\n")
  cat(paste("The estimation took:", make.time(x) ,"\n"))
  
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits)
  
  cat(paste("\nOptimization of log-likelihood by", x$logLik$type))
  cat(paste("\nLog Likelihood:", signif(x$logLik$maximum, digits)))
  cat(paste("\nNumber of observations:",  x$logLik$nobs))
  cat(paste("\nNumber of iterations:" , x$logLik$iterations))
  cat(paste("\nExit of MLE:", x$logLik$message))
  
  if (!(x$model == "mnl" | x$model == "lc")) cat(paste("\nSimulation based on", x$R, "draws"))
  invisible(x)
}

#' vcov method for gmnl objects
#' 
#' The \code{vcov} method for \code{gmnl} objects extracts the covariance matrix of the coefficients or the random parameters. It also allows to get the standard errors for the variance-covariance matrix of the random parameters
#' 
#' @param object a fitted model of class \code{gmnl},
#' @param what indicates which covariance matrix has to be extracted. The default is \code{coefficient}, in this case the \code{vcov} behaves as usual. If \code{what = "ranp"} the covariance matrix of the random parameters is returned as default, 
#' @param type if the model is estimated with random parameters, then this argument indicates what matrix should be returned. If \code{type = "cov"}, then the covariance matrix of the random parameters is returned; if \code{type = "cor"} then the correlation matrix of the random parameters is returned; if \code{type = "sd"} then the standard deviation of the random parameters is returned,
#' @param se if \code{TRUE} \code{type = "cov"} then the standard error of the covariance matrix of the random parameters is returned; if \code{TRUE} \code{type = "sd"} the standard error of the standard deviation of the random parameter is returned. This argument if valid only if the model is estimated using correlated random parameters,
#' @param Q this argument is only valid if the "\code{mm}" (MM-MNL) model is estimated. It indicates the class for which the variance-covariance matrix is computed,
#' @param digits number of digits,
#' @param ... further arguments
#' @details This new interface replaces the \code{cor.gmnl}, \code{cov.gmnl} and \code{se.cov.gmnl} functions which are deprecated.
#' @seealso \code{\link[gmnl]{gmnl}} for the estimation of multinomial logit models with random parameters.
#' @method vcov gmnl
#' @import stats
#' @export
vcov.gmnl <- function(object, what = c('coefficient', 'ranp'), type = c('cov', 'cor', 'sd'), 
                      se =  FALSE, Q = NULL, digits = max(3, getOption("digits") - 2), ...)
{
  what <- match.arg(what)
  type <- match.arg(type)
  if (what == 'coefficient') {
    H    <- object$logLik$hessian
    vcov <- solve(-H)
    rownames(vcov) <- colnames(vcov) <- names(coef(object))
    return(vcov)
  }
  if (what == 'ranp') {
    if (se) {
      if (type == 'cov') se.cov.gmnl(object, sd = FALSE, Q = Q, digits = digits)
      if (type == 'sd')  se.cov.gmnl(object, sd = TRUE, Q = Q,  digits = digits)
      if (type == 'cor') stop("standard error for correlation coefficients not implemented yet")
    } else {
      if (type == 'cov') print(cov.gmnl(object, Q = Q)) 
      if (type == 'cor') print(cor.gmnl(object, Q = Q))
      if (type == 'sd')  print(sqrt(diag(cov.gmnl(object, Q))))
    }

  }
}

#' @rdname gmnl
#' @method update gmnl
#' @import stats
#' @export
update.gmnl <- function(object, new, ...){
  call <- object$call
  if (is.null(call))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(new))
    call$formula <- update(formula(object), new)
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, parent.frame())
}

#' @rdname gmnl
#' @export 
coef.gmnl <- function(object, ...){
  result <- object$coefficients
  return(result)
}

#' @rdname gmnl
#' @export
model.matrix.gmnl <- function(object, ...){
  model.matrix(object$formula, object$mf)
}


model.response.gmnl <- function(object, ...){
  y.name <- paste(deparse(object$formula[[2]]))
  object$mf[[y.name]]
}

#' @rdname gmnl
#' @export
residuals.gmnl <- function(object, outcome = TRUE, ...){
  if (!outcome) {
    result <- object$residuals
  }
  else{
    J <- ncol(object$residuals)
    y <- matrix(model.response.gmnl(object), ncol = J, byrow = T)
    result <- apply(y * object$residuals, 1, sum)
  }
  result
}

#' @rdname gmnl
#' @import stats
df.residual.gmnl <- function(object, ...){
  n <- length(residuals(object))
  K <- length(coef(object))
  return(n - K)
}

#' @rdname gmnl
#' @export
fitted.gmnl <- function(object, outcome = TRUE, ...){
  if (outcome) result <- object$prob.ind
  else result <- object$prob.alt
  result
}

#' @rdname gmnl
#' @export
logLik.gmnl <- function(object,...){
  structure(object$logLik$maximum[[1]], df = length(object$coefficients),
            nobs = object$logLik$nobs, class = "logLik")
}


#' Get Model Summaries for Use with "mtable"
#' 
#' A generic function to collect coefficients and summary statistics from a \code{gmnl} object. It is used in \code{mtable}.
#' 
#' @param obj a \code{gmnl} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}
#' @examples
#' ## Estimate MNL models
#' data("TravelMode", package = "AER")
#' library(mlogit)
#' TM <- mlogit.data(TravelMode, choice = "choice", shape = "long", 
#'                  alt.levels = c("air", "train", "bus", "car"), chid.var = "individual")
#'                  
#' mnl.1 <- gmnl(choice ~ wait + vcost + travel + gcost | 0, data = TM)
#' mnl.2 <- gmnl(choice ~ wait + vcost                  | 0, data = TM) 
#' 
#' ## Table
#' library(memisc)
#' mtable("MNL 1"= mnl.1, "MNL 2" = mnl.2, 
#'        summary.stats = c("N", "Log-likelihood", "BIC", "AIC"))
#' @import stats
#' @export getSummary.gmnl
getSummary.gmnl <- function(obj, alpha = 0.05, ...){
  smry <- summary(obj)
  coef <- smry$CoefTable
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha / 2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha / 2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  obj$logLik$nobs
  ll <- logLik(obj)
  sumstat <- c(logLik = ll, deviance = NA, AIC = AIC(obj), BIC = BIC(obj), N = N, 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
               Nagelkerke = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}

#' Akaike's Information Criterion
#' 
#' Calculate the Akaike's information Criterion (AIC) or the Bayesian
#' information Criterion (BIC) for an object of class \code{gmnl}.
#' 
#' @param object a fitted model of class \code{gmnl}.
#' @param ... additional arguments to be passed to or from other functions.
#' @param k a numeric value, use as penalty coefficient for number of parameters
#' in the fitted model.
#' @details For more information see \code{\link[stats]{AIC}} or \code{\link[stats]{BIC}}
#' @return A numeric value with the corresponding AIC or BIC value.
#' @seealso \code{\link[gmnl]{gmnl}} for the estimation of multinomial logit models with observed and unobserved individual heterogeneity.
#' 
#' @import stats
#' @method AIC gmnl
#' @export
AIC.gmnl <- function(object, ..., k = 2){
  return(-2 * object$logLik$maximum[[1]] + k * length(coef(object)))
}

#' @rdname AIC.gmnl
#' @import stats
#' @method BIC gmnl
#' @export 
#' @examples
#' 
#' ## Estimate MNL model
#' data("TravelMode", package = "AER")
#' library(mlogit)
#' TM <- mlogit.data(TravelMode, choice = "choice", shape = "long", 
#'                  alt.levels = c("air", "train", "bus", "car"), chid.var = "individual")
#'                  
#' mnl <- gmnl(choice ~ wait + vcost + travel + gcost | 0 , data = TM)
#' AIC(mnl)
#' BIC(mnl)
BIC.gmnl <- function(object, ...){
  return(AIC(object, k = log(object$logLik$nobs)))
}

#### Methods for sandwiches

#' Bread for Sandwiches
#' 
#' Computes the ``bread'' of the sandwich covariance matrix for objects of class \code{gmnl}.
#' 
#' @param x a fitted model of class \code{gmnl}.
#' @param ... other arguments when \code{bread} is applied to another
#' class object.
#' @return The covariance matrix times observations
#' @details For more information see \code{\link[sandwich]{bread}} from the package \pkg{sandwich}.
#' @references Zeileis A (2006), Object-oriented Computation of Sandwich 
#' Estimators. Journal of Statistical Software, 16(9), 1--16.
#' @method bread gmnl
#' @import stats
#' @export bread.gmnl
bread.gmnl <- function(x, ... ){
  return( vcov( x ) * x$logLik$nobs)
}

#' Gradient for Observations
#' 
#' It extracts the gradient for each observation evaluated at the estimated parameters for an object of class \code{gmnl}.
#' 
#' @param x a fitted model of class \code{gmnl}.
#' @param ... other arguments. Ignored.
#' @return The gradient matrix of dimension \eqn{n \times K}
#' @references Zeileis A (2006), Object-oriented Computation of Sandwich 
#' Estimators. Journal of Statistical Software, 16(9), 1--16.
#' @details For more information see \code{\link[sandwich]{estfun}} from package \pkg{sandwich}.
#' @method estfun gmnl
#' @export estfun.gmnl
estfun.gmnl <- function(x, ... ){
  return(x$logLik$gradientObs )
}

#' @rdname gmnl
#' @export nObs.gmnl
nObs.gmnl <- function(x, ... ){
  return(x$logLik$nobs)
}

#' Get the Conditional Individual Coefficients
#' 
#' This a helper function to obtain the individuals' conditional estimate of the either random parameters or willingness-to-pay.
#' @param x an object of class \code{gmnl}.
#' @param par a string giving the name of the variable with a random parameter.
#' @param effect a string indicating what should be computed: the conditional expectation of the individual coefficients "\code{ce}", or the conditional expectation of the willingness-to-pay "\code{wtp}".
#' @param wrt a string indicating with respect to which variable the willingness-to-pay should be computed.
#' @param ... further arguments. Ignorred.
#' 
#' @return A named list where "\code{mean}" contains the individuals' conditional mean for the random parameter or willingness-to-pay, and where "\code{sd.est}" contains standard errors. 
#' @export
#' @author Mauricio Sarrias.
#' @references
#' \itemize{
#' \item Greene, W. H. (2012). Econometric Analysis, Seventh Edition. Pearson Hall.
#' \item Train, K. (2009). Discrete Choice Methods with Simulation. Cambridge University Press.
#' }
#' @seealso \code{\link[gmnl]{gmnl}} for the estimation of multinomial Logit models with individual parameters.
#' @import stats
#' @examples
#' \dontrun{
#' ## Data
#' data("TravelMode", package = "AER")
#' library(mlogit)
#' TM <- mlogit.data(TravelMode, choice = "choice", shape = "long", 
#'                  alt.levels = c("air", "train", "bus", "car"), chid.var = "individual")
#'                  
#' ## MIXL model with observed heterogeneity
#' mixl.hier <- gmnl(choice ~ vcost + gcost + travel + wait | 1 | 0 | income + size - 1,
#'                  data = TM,
#'                  model = "mixl",
#'                  ranp = c(travel = "t", wait = "n"),
#'                  mvar = list(travel = c("income","size"), wait = c("income")),
#'                  R = 30,
#'                  haltons = list("primes"= c(2, 17), "drop" = rep(19, 2)))
#'                  
#' ## Get the individuals' conditional mean and their standard errors for lwage
#' bi.travel <- effect.gmnl(mixl.hier, par = "travel", effect = "ce")
#' summary(bi.travel$mean)
#' summary(bi.travel$sd.est)
#' 
#' ## Get the individuals' conditional WTP of travel with respect to gcost
#' wtp.travel <- effect.gmnl(mixl.hier, par = "travel", effect = "wtp", wrt = "gcost")
#' summary(wtp.travel$mean)
#' summary(wtp.travel$sd.est)
#' } 
effect.gmnl <- function(x, par = NULL, effect = c("ce", "wtp"), wrt = NULL, ... ){
  if (!inherits(x, "gmnl")) stop("not a \"gmnl\" object")
  model <- x$model
  if (model == "mnl") stop("This function is valid only for models with individual heterogeneity")
  type <- match.arg(effect)
  ranp <- x$ranp
  #if (model != "lc" && !is.null(par) && !(par %in% names(ranp))) stop("This parameter is not random: ", par)
  #if (model != "lc" ||  model!= "smnl") if (!(par %in% names(ranp))) stop("This parameter is not random: ", par)
  if (type == "wtp" & is.null(wrt)) stop("you need to specify wrt")
  bi   <- x$bi
  Qir <- x$Qir
  if (model == "mixl" || model == "gmnl" || model == "smnl") {
    N <- nrow(Qir)
    K <- dim(bi)[[3]]
    var_coefn <- dimnames(bi)[[3]]
    mean <- mean.sq <- matrix(NA, N, K)
    if (type == "wtp") {
      if (model != "smnl") {
        is.ran <- any(names(ranp) %in% wrt)
        gamma <- if (is.ran) bi[, , wrt] else coef(x)[wrt]
      } else gamma <- bi[, , wrt]
      for (j in 1:K) {
        mean[, j]    <- rowSums((bi[, , j] / gamma)   * Qir)
        mean.sq[, j] <- rowSums(((bi[, , j] / gamma) ^ 2) * Qir)
      }
    } else {
      for (j in 1:K) {
        mean[, j]    <-  rowSums(bi[, , j] * Qir)
        mean.sq[, j] <-  rowSums(bi[, , j] ^ 2 * Qir) 
      }
    }
  }
  if (model == "lc") {
    N <- nrow(Qir)
    K <- ncol(bi)
    var_coefn <- colnames(bi)
    mean <- mean.sq <- matrix(NA, N, K)
    if (type == "wtp") {
      gamma <- bi[, wrt]
      for (j in 1:K) {
        mean[, j]    <- rowSums(repRows(bi[, j] / gamma, N) * Qir)
        mean.sq[, j] <- rowSums(repRows((bi[, j] / gamma) ^ 2, N) * Qir)
      }
    } else {
      for (j in 1:K) {
        mean[, j]    <- rowSums(repRows(bi[, j], N)   * Qir) 
        mean.sq[, j] <- rowSums(repRows(bi[, j] ^ 2, N) * Qir) 
      }
    }
  }
  if (model == "mm") {
    wnq  <- Qir$wnq
    Ln   <- Qir$Ln
    Pnrq <- Qir$Pnrq
    N <- length(Ln)
    K <- dim(bi)[[4]]
    mean <- mean.sq <- matrix(NA, N, K)
    var_coefn <- dimnames(bi)[[4]]
    if (type == "wtp") {
      gamma <- bi[,,,wrt]
      for (j in 1:K) {
        mean[, j]    <- rowSums(wnq * apply((bi[,,,j] / gamma)   * Pnrq, c(1, 3), mean) / Ln)
        mean.sq[, j] <- rowSums(wnq * apply((bi[,,,j] / gamma) ^ 2 * Pnrq, c(1, 3), mean) / Ln)
      }
    } else {
      for (j in 1:K) {
        mean[, j]    <- rowSums(wnq * apply(bi[,,,j]   * Pnrq, c(1, 3), mean) / Ln)
        mean.sq[, j] <- rowSums(wnq * apply(bi[,,,j] ^ 2 * Pnrq, c(1, 3), mean) / Ln)
      }
    }
  }
  sd.est  <- suppressWarnings(sqrt(mean.sq - mean ^ 2))
  colnames(mean) <- colnames(sd.est) <- var_coefn
  if (!is.null(par)) {
    mean   <- mean[, par]
    sd.est <- sd.est[, par]
  }
  effe <- list(
               mean = mean,
               sd.est = sd.est)
  return(effe)
}


#' Plot of the Distribution of the Conditional Expectation of Random Parameters
#' 
#' Methods for \code{gmnl} objects which provide a plot of the distribution of the conditional expectation of the random parameters or the distribution of the conditional willigness-to-pay.
#' 
#' 
#' @param x an object of class \code{gmnl}.
#' @param par a string giving the name of the variable with random parameter.
#' @param type a string indicating the type of distribution: it can be a \code{histogram} or a \code{density} of the conditional expectation of the random coefficients or WTP.
#' @param ind a boolean. If \code{TRUE}, a 95\% interval of conditional distribution for each individual is plotted. As default, the conditional expectation of \code{par} for the first 10 individual is plotted.
#' @param id only relevant if \code{ind} is not \code{NULL}. This is a vector indicating the individuals for whom the user want to plot the conditional coefficients.
#' @param effect a string indicating whether the conditional expectation, "\code{ce}", or the WTP, "\code{wtp}" should be plotted.
#' @param wrt a string indicating with respect to which variable the WTP should be computed if \code{effect = "wtp"}. 
#' @param adjust  bandwidth for the kernel density.
#' @param main an overall title for the plot.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param col color for the graph.
#' @param breaks number of breaks for the histrogram if \code{type = "histogram"}.
#' @param ... further arguments to be passed to \code{plot} or \code{plotCI}.
#' @references
#' \itemize{
#' \item Greene, W. H. (2012). Econometric Analysis, Seventh Edition. Pearson Hall.
#' \item Train, K. (2009). Discrete Choice Methods with Simulation. Cambridge University Press.
#' }
#' @seealso \code{\link[gmnl]{gmnl}} for the estimation of different multinomial models with individual heterogeneity and \code{\link[gmnl]{effect.gmnl}}.
#' @importFrom plotrix plotCI
#' @method plot gmnl
#' @author Mauricio Sarrias
#' @export 
#' @import graphics
#' @import stats
#' @examples
#' \dontrun{
#' ## Examples using the Electricity data set from the mlogit package
#' library(mlogit)
#' data("Electricity", package = "mlogit")
#' Electr <- mlogit.data(Electricity, id.var = "id", choice = "choice",
#'                      varying = 3:26, shape = "wide", sep = "")
#'                      
#' ## Estimate a MIXL model with correlated random parameters
#' Elec.cor <- gmnl(choice ~ pf + cl + loc + wk + tod + seas| 0, data = Electr,
#'                  subset = 1:3000,
#'                  model = 'mixl',
#'                  R = 10,
#'                  panel = TRUE,
#'                  ranp = c(cl = "n", loc = "n", wk = "n", tod = "n", seas = "n"),
#'                  correlation = TRUE)
#'                  
#' ## Plot the density of the conditional expectation distribution of loc
#' plot(Elec.cor, par = "loc", effect = "ce", type = "density", col = "grey")
#' 
#' ## Plot the conditional expectation of loc for each individual
#' plot(Elec.cor, par = "loc", effect = "ce", ind = TRUE, id = 1:30)
#' 
#' ## Plot the WTP for cl
#' plot(Elec.cor, par = "loc", effect = "wtp", wrt = "pf")                  
#'} 
plot.gmnl <- function(x, par = NULL, effect = c("ce", "wtp"), wrt = NULL,
                              type = c("density", "histogram"), adjust = 1, 
                              main = NULL, col = "indianred1", breaks = 10, ylab = NULL,
                              xlab = NULL, ind = FALSE, id = NULL, ...){
  model <- x$model
  if (model == "mnl") stop("The plot is valid only for models with individual heterogeneity")
  if (is.null(par)) stop("Must specified the name of the parameter")
  type   <- match.arg(type)
  effect <- match.arg(effect)
  xlab <- switch(effect,
                 "wtp" = expression(E(hat(wtp[i]))),
                 "ce"  = expression(E(hat(beta[i])))) 
  if (!ind) {
    if (is.null(main)) main <- paste("Conditional Distribution for", par)
    if (is.null(ylab)) {
      ylab <- switch(type,
                     "density"   = "Density",
                     "histogram" = "Frequency")
    }
    rpar <- effect.gmnl(x, par,  effect = effect, wrt =  wrt)$mean
    if (type == "density") {
      pdens <- density(rpar, adjust = adjust)
      plot(pdens, ylab = ylab, xlab = xlab, main = main, col = col)
      has.pos <- any(pdens$x > 0)
      if (has.pos) {
        x1 <- min(which(pdens$x >= 0))  
        x2 <- max(which(pdens$x <  max(pdens$x)))
        with(pdens, polygon(x = c(x[c(x1, x1:x2, x2)]), y = c(0, y[x1:x2], 0), col = col, border = NA))
      }
    } else {
      minb <- round(min(rpar), 2)
      maxb <- round(max(rpar), 2)
      hist(rpar, xlab = xlab, main = main, col = col, breaks = breaks, 
           xaxs = "i", yaxs = "i", las = 1, xaxt = 'n', ylab = ylab)
      axis(1, at = seq(minb, maxb, (maxb - minb) * .05))
    }
  } else {
    if (is.null(main)) main <- paste("95% Probability Intervals for ", par)
    if (is.null(id)) id <- seq(1, 10, 1) 
    if (is.null(ylab)) ylab <- "Individuals"
    f.bran <- effect.gmnl(x, par,  effect = effect, wrt =  wrt)$mean
    f.sran <- effect.gmnl(x, par,  effect = effect, wrt =  wrt)$sd.est
    lower <- f.bran - qnorm(0.975) * f.sran
    upper <- f.bran + qnorm(0.975) * f.sran
    plotrix::plotCI(as.numeric(id), f.bran[id], ui = upper[id], li = lower[id],
                    xlab = ylab, ylab = xlab,
                    lty = 2, main = main,
                    pch = 21, col = col)
  } 
}



#' Functions for Correlated Random Parameters
#' 
#' These are a set of functions that help to extract the variance-covariance matrix, the correlation matrix, and the standard error of the random parameters for models of class \code{gmnl}.
#' 
#' @param x an object of class \code{gmnl} where \code{ranp} is not \code{NULL}.
#' @param Q this argument is only valid if the "\code{mm}" (MM-MNL) model is estimated. It indicates the class for which the variance-covariance matrix is computed.
#' @param sd if \code{TRUE}, then the standard deviations of the random parameters along with their standard errors are computed. 
#' @param digits the number of digits.
#' 
#' @return \code{cov.gmnl} returns a matrix with the variance of the random parameters if the model is fitted with random coefficients. If the model is fitted with \code{correlation = TRUE}, then the variance-covariance matrix is returned. 
#' 
#'   
#' If \code{correlation = TRUE} in the fitted model, then  \code{se.cov.gmnl} returns a coefficient matrix for the elements of the variance-covariance matrix or the standard deviations if \code{sd = TRUE}.
#' 
#' 
#' @details The variance-covariance matrix is computed using the Cholesky decomposition \eqn{LL'=\Sigma}.
#' 
#' 
#' \code{se.cov.gmnl} function is a wrapper for the \code{\link[msm]{deltamethod}} function of the \pkg{msm} package.
#' @author Mauricio Sarrias \email{msarrias86@@gmail.com}
#' @references
#' \itemize{
#' \item Greene, W. H. (2012). Econometric Analysis, Seventh Edition. Pearson Hall.
#' \item Train, K. (2009). Discrete Choice Methods with Simulation. Cambridge University Press.
#' }
#' @seealso \code{\link[gmnl]{gmnl}} for the estimation of different multinomial models with individual heterogeneity.
#' @examples
#' \dontrun{
#' ## Examples using Electricity data set from mlogit package
#' library(mlogit)
#' data("Electricity", package = "mlogit")
#' Electr <- mlogit.data(Electricity, id.var = "id", choice = "choice",
#'                      varying = 3:26, shape = "wide", sep = "")
#'                      
#' ## Estimate a MIXL model with correlated random parameters
#' Elec.cor <- gmnl(choice ~ pf + cl + loc + wk + tod + seas| 0, data = Electr,
#'                  subset = 1:3000,
#'                  model = 'mixl',
#'                  R = 10,
#'                  panel = TRUE,
#'                  ranp = c(cl = "n", loc = "n", wk = "n", tod = "n", seas = "n"),
#'                  correlation = TRUE)
#'                  
#' ## Use functions for correlated random parameters
#' cov.gmnl(Elec.cor)
#' se.cov.gmnl(Elec.cor)
#' se.cov.gmnl(Elec.cor, sd = TRUE)
#' cor.gmnl(Elec.cor)
#' }
#' @export
cov.gmnl <- function(x, Q = NULL){
  if (!inherits(x, "gmnl")) stop("not a \"gmnl\" object")
  if (is.null(x$ranp)) stop('cov.gmnl only relevant for random coefficient model')
  model <- x$model
  if (!is.null(Q) & model != "mm") stop("Q is only relevant for MM-MNL model")
  if (model == "mm") {
    if (is.null(Q)) stop("MM-MNL model requires Q")
    if (Q > x$Q) stop("Q is greater than the number of classes in the fitted model")
  }  
  beta.hat <- x$coefficients
  K  <- length(x$ranp)
  nr <- names(x$ranp)
  if (x$correlation) {
    names.stds <- c()
    if (model == "mm") {
      for (i in 1:K) names.stds <- c(names.stds, paste('class', Q, 'sd', nr[i], nr[i:K], sep = '.'))
    } else {
      for (i in 1:K) names.stds <- c(names.stds, paste('sd', nr[i], nr[i:K], sep = '.'))
    }
    v    <- beta.hat[names.stds]
    V    <- tcrossprod(makeL(v))
    colnames(V) <- rownames(V) <- nr
  } else{
    names.stds <- if (model != "mm") paste("sd", nr, sep = ".") else paste("class", Q, "sd", nr, sep = ".")
    sv   <- beta.hat[names.stds]
    V    <- matrix(0, K, K)
    diag(V) <- sv ^ 2
    colnames(V) <- rownames(V) <- nr
  }
  V
}

#' @rdname cov.gmnl
#' @export
cor.gmnl <- function(x, Q = NULL){
  if (!x$correlation) stop('cor.gmnl only relevant for correlated random coefficient')
  V   <- cov.gmnl(x, Q = Q)
  nr  <- names(x$ranp)
  D   <- diag(sqrt(diag(V)))
  Rho <- solve(D) %*% V %*% solve(D)
  colnames(Rho) <- rownames(Rho) <- nr
  Rho
}


#' @rdname cov.gmnl
#' @importFrom msm deltamethod
#' @import stats
#' @export
se.cov.gmnl <- function(x, sd =  FALSE, Q = NULL, digits = max(3, getOption("digits") - 2)){
  if (!inherits(x, "gmnl")) stop("not a \"gmnl\" object")
  if (!x$correlation) stop('se.cov.gmnl only relevant for correlated random coefficient')
  model <- x$model
  if (!is.null(Q) & model != "mm") stop("Q is only relevant for MM-MNL model")
  if (model == "mm") {
    if (is.null(Q)) stop("MM-MNL model requires Q")
    if (Q > x$Q) stop("Q is greater than the number of classes in the fitted model")
  }  
  beta.hat <- x$coefficients
  Ka <- length(x$ranp)
  nr <- names(x$ranp)
  names.stds <- c()
  if (model == "mm") {
    for (i in 1:Ka) names.stds <- c(names.stds, paste('class', Q, 'sd', nr[i], nr[i:Ka], sep = '.'))
  } else {
    for (i in 1:Ka) names.stds <- c(names.stds, paste('sd', nr[i], nr[i:Ka], sep = '.'))
  }
  stds.hat <- beta.hat[names.stds]
  sel.vcov <- vcov(x)[names.stds, names.stds]
  form <- c()
  if (sd) {
    for (i in 1:Ka) {
      k <- i
      if (i == 1) {
        form <- paste("~ sqrt(", c(form, paste(paste("x",  i, sep = ""), paste("x", k, sep = ""), sep = "*")), ")")
      } else {
        temp <- paste(paste("x",  i, sep = ""), paste("x", k, sep = ""), sep = "*")
        j <- 2
        while(j <= i) {
          temp <- paste(temp, make.add(row = j, col = k, Ka = Ka)[1], sep = "+") 
          j <- j + 1
        }
        form <- c(form, paste("~ sqrt(", temp, ")"))
      }
    }
    b <- sqrt(diag(cov.gmnl(x, Q)))
    names(b) <- colnames(cov.gmnl(x, Q))
  } else {
    for (i in 1:Ka) { 
      if (i == 1) {
        form <- paste("~", c(form, paste(paste("x",  i:Ka, sep = ""), paste("x", i, sep = ""), sep = "*")))
      } else {
        temp <- paste(paste("x",  i:Ka, sep = ""), paste("x", i, sep = ""), sep = "*")
        j <- 2
        while(j <= i) {
          temp <- paste(temp, make.add(row = j, col = i, Ka = Ka), sep = "+") 
          j <- j + 1
        }
        form <- c(form, paste("~", temp))
      }
    }
    names.vcov <- c()
    for (i in 1:Ka) names.vcov <- c(names.vcov, paste('v', nr[i], nr[i:Ka], sep = '.'))
    b <- drop(cov.gmnl(x, Q)[lower.tri(cov.gmnl(x, Q), diag = TRUE)])
    names(b) <- names.vcov
  }
  std.err <- c()
  for (i in 1:length(form)) {
    std.err <- c(std.err, msm::deltamethod(as.formula(form[i]), stds.hat, sel.vcov, ses =  TRUE))
  }  
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  tableChol <- cbind(b, std.err, z, p)
  if (!sd) cat(paste("\nElements of the variance-covariance matrix \n\n"))
  else cat(paste("\nStandard deviations of the random parameters \n\n"))
  #colnames(tableChol) <- c("Estimate", "Std. Error", "t-value", "Pr(>|t|)") 
  colnames(tableChol) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)") 
  printCoefmat(tableChol, digits =  digits)
}

#' Compute Willingness-to-pay
#' 
#' Compute the willingness-to-pay.
#' 
#' @param object an object of class \code{gmnl}.
#' @param wrt a string indicating the variable with respect to which the WTP is computed,
#' @param digits number of significant digits to be used for most numbers.
#' @return A coefficient matrix with the WTP point estimates and standard errors. 
#' @export
#' @details For each coefficient, this function computes both the point estimate and standard error of WTP with respect to the variable specified in the argument \code{wrt}. Specifically, let \eqn{\beta_k} be the coefficient for variable \eqn{k}, then \deqn{WTP_{k}=-\beta_k/\beta_p}
#' 
#' 
#' where \eqn{\beta_p} is the coefficient for the variable specified with the argument \code{wrt}. Note that, \code{wtp.gmnl} does not include the negative sign. 
#' 
#' 
#' \code{wtp.gmnl} function is a wrapper for the \code{\link[msm]{deltamethod}} function of the \pkg{msm} package. 
#' @seealso \code{\link[msm]{deltamethod}} for the estimation of the standard errors.
#' @author Mauricio Sarrias.
#' @examples
#' 
#' ## Examples using the Electricity data set from the mlogit package
#' library(mlogit)
#' data("Electricity", package = "mlogit")
#' Electr <- mlogit.data(Electricity, id.var = "id", choice = "choice",
#'                      varying = 3:26, shape = "wide", sep = "")
#'                      
#' ## Estimate a conditional logit model
#' clogit <- gmnl(choice ~ pf + cl + loc + wk + tod + seas| 0,
#'                data = Electr)
#' wtp.gmnl(clogit, wrt = "pf")
#' @import stats
#' @references
#' \itemize{
#' \item Greene, W. H. (2012). Econometric Analysis, Seventh Edition. Pearson Hall.
#' \item Train, K. (2009). Discrete Choice Methods with Simulation. Cambridge University Press.
#' }
wtp.gmnl <- function(object, wrt =  NULL, digits = max(3, getOption("digits") - 2)){
  if (is.null(wrt)) stop("WTP needs the variable in the denominator: wrt")
  beta.hat <- coef(object)
  posi <- match(wrt, names(beta.hat))
  form <- c()
  b <- c()
  namesb <- names(beta.hat)[-c(posi)]
  for (i in 1:length(beta.hat)) {
    if (i != posi) {
      b <- c(b, beta.hat[i]/ beta.hat[posi])
      form <- c(form, paste("~", "x", i, "/", "x", posi, sep = ""))
    }
  }
  names(b) <- namesb
  std.err <- c()
  for (i in 1:length(form)) {
    std.err <- c(std.err, msm::deltamethod(as.formula(form[i]), beta.hat, vcov(object), ses =  TRUE))
  }
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  tablewtp <- cbind(b, std.err, z, p)
  colnames(tablewtp) <- c("Estimate", "Std. Error", "t-value", "Pr(>|t|)") 
  cat(paste("\nWilligness-to-pay respect to: ", wrt, "\n\n"))
  printCoefmat(tablewtp, digits = digits)
}
