#' Fitting GWRM Models
#'
#' \code{gw} is used to fit Generalized Waring Regression Models (GWRM), specified by giving a symbolic description of the linear predictor.
#' 
#' @aliases gw gw.fit
#'
#' @param formula	an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data	an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}.
#' @param weights	an optional vector of 'prior weights' to be used in the fitting process. Should be \code{NULL} or a numeric vector.
#' @param k	optional value for the \code{k} parameter. If \code{NULL}, it is estimated.
#' @param subset	an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action	a function which indicates what should happen when the data contain \code{NA} values. See \code{\link{glm}}.
#' @param kstart	starting value for the \code{k} parameter.
#' @param rostart	starting value for the \code{ro} parameter.
#' @param betastart	starting values for the vector of means.
#' @param offset	this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be \code{NULL} or a numeric vector of length equal to the number of cases. One or more offset terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.
#' @param control	a list of parameters for controlling the fitting process.
#' @param method	the method to be used in fitting the model. The default method initially uses non-linear minimization (\code{nlm}) and Nelder-Mead optimization (\code{optim}) to fit a model which is then re-fitted by \code{"L-BFGS-B"} (\code{optim}). In this way, SE estimates for all the model parameters are provided. \code{"nlm"} and \code{"Nelder-Mead"} are also possible values, but they do not provide SE estimates for \code{k} and \code{ro}.
#' @param hessian	if \code{TRUE}, the hessian of \code{f} at the minimum is returned.
#' @param model	a logical value indicating whether model frame should be included as a component of the returned value.
#' @param x,y  For \code{gw}: logical values indicating whether the response vector and model matrix used in the fitting process should be returned as components of the returned value.
#' 
#' For \code{gw.fit}: \code{x} is a design matrix of dimension \code{n * p}, and \code{y} is a vector of observations of length \code{n}.
#' @param ...	further arguments.
#'
#' @return \code{gw} returns an object of class \code{"gw"}. The function \code{summary} can be used to obtain or print a summary of the results. An object of class \code{"gw"} is a list containing the following components:
#' \itemize{
#' \item \code{Y} {if requested (the default), the \code{y} vector used.}
#' \item \code{W} {the weights supplied, a vector of \code{1}s if none were.}
#' \item \code{covars} {names of the covariates in the model.}
#' \item \code{nobs} {number of observations.}
#' \item \code{covoffset} {a logical value specifying if an offset is present.}
#' \item \code{loglik} {the maximized log-likelihood.}
#' \item \code{aic} {a version of Akaike's \emph{An Information Criterion}, minus twice the maximized log-likelihood plus twice the number of parameters.}
#' \item \code{bic} {Bayesian Information Criterion, minus twice the maximized log-likelihood plus the number of parameters multiplied by the logarithm of the number of observations.}
#' \item \code{df.residual} {the residual degrees of freedom.}
#' \item \code{residuals} {the residuals in the final iteration of the fit.}
#' \item \code{coefficients} {a named vector of coefficients.}
#' \item \code{betaIIpars} {parameters estimates of the BetaII distribution.}
#' \item \code{betascoefs} {a vector of coefficients.}
#' \item \code{fitted.values} {the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
#' \item \code{hessian} {a symmetric matrix giving an estimate of the Hessian at the solution found in the optimization of the log-likelihood function.}
#' \item \code{cov} {an estimate of the covariance matrix of the model coefficients.}
#' \item \code{se} {a vector of the standard errors estimates of the estimated coefficients.}
#' \item \code{corr} {an estimate of the correlation matrix of the model coefficients.}
#' \item \code{code} {a code that indicates successful convergence of the fitter function used (see \code{nlm} and \code{optim} helps).}
#' \item \code{converged} {logical value that indicates if the optimization algorithms succesfull.}
#' \item \code{method} {the name of the fitter function used.}
#' \item \code{k} {if requested, the \code{k} value used.}
#' \item \code{kBool} {a logical value specifying whether there is a \code{k} value or it is estimated.}
#' \item \code{call} {the matched call.}
#' \item \code{formula} {the formula supplied.}
#' \item \code{terms} {the \code{terms} object used.}
#' \item \code{data} {the \code{data} argument.}
#' \item \code{offset} {the offset vector used.}
#' \item \code{control} {the value of the \code{control} argument used.}
#' \item \code{method} {the name of the fitter function used.}
#' \item \code{contrasts} {(where relevant) the contrasts used.}
#' \item \code{xlevels} {(where relevant) a record of the levels of the factors used in fitting.}
#' }
#'
#' @importFrom stats model.response is.empty.model model.matrix contrasts model.weights model.offset AIC .getXlevels
#'
#' @examples
#' data(goals)
#' gw(goals ~ position + offset(log(played)), data = goals)
#'
#' @export

gw <- function(formula, data, weights, k = NULL, subset, na.action,
               kstart = 1, rostart = 2, betastart = NULL,offset, control = list(...), method = NULL, hessian = TRUE, model = TRUE, x = FALSE, y = TRUE, ...){

  warningDefault<-getOption("warn")

  if (is.null(control$trace)){
    control$trace = 0
  }

  if (control$trace== 0)
    options(warn=-1)

  call <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action","offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval.parent(mf)

  fitted<-FALSE

  Terms <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  nobs<-nrow(as.matrix(Y))
  X <- if(!is.empty.model(Terms)){
    model.matrix(Terms, mf, contrasts)
  }
  else{
    matrix(1, nobs, 0)
  }
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))

  if (!is.null(offset)) {
    if (length(offset) != nobs)
      stop(gettextf("Number of offsets is %d should equal %d (number of observations)", length(offset), nobs), domain = NA)
  }
  if (is.null(method)){
    if (control$trace > 0) {
      cat("Trying 'nlm' initial fit", "\n")
    }
    nlm.fit <-try(eval(gw.fit(x = X, y = Y, weights = weights, k = k, kstart = kstart, rostart = rostart, betastart = betastart, offset = offset, control = control, method = "nlm", hessian = hessian)), silent = TRUE)
    if ('try-error' %in% class(nlm.fit)){
      if (control$trace > 0) {
        cat("Crashed 'nlm' initial fit", "\n")
      }
      nlm.aic <- Inf
    }
    else {
      nlm.aic <- AIC(nlm.fit)
      if (control$trace > 0) {
        cat(paste(round(nlm.aic, 2), "AIC value in 'nlm' initial fit"), "\n")
      }
    }
    if (control$trace > 0) {
      cat("Trying 'Nelder-Mead' initial fit", "\n")
    }
    neldermead.fit <- try(eval(gw.fit(x = X, y = Y, weights = weights, k = k, kstart = kstart, rostart = rostart, betastart = betastart, offset = offset, control = control, method = "Nelder-Mead", hessian = hessian)), silent = TRUE)

    if ('try-error' %in% class(neldermead.fit)){
      if (control$trace > 0) {
        cat("Crashed 'Nelder-Mead' initial fit", "\n")
      }
      neldermead.aic <- Inf
    }
    else {
      neldermead.aic <- AIC(neldermead.fit)
      if (control$trace > 0) {
        cat(paste(round(neldermead.aic,2), "AIC value in 'Nelder-Mead' initial fit"), "\n")
      }
    }
    if (neldermead.aic == Inf && nlm.aic == Inf){
      warning("No 'nlm' neither 'Nelder-Mead' provide fits. Try to change initial values")
    }
    else {
      if (neldermead.aic < nlm.aic | nlm.aic <= 0){
        fit2 <- neldermead.fit
      }
      else{
        fit2<-nlm.fit
      }
      if (control$trace > 0) {
        cat("L-BFGS-B fitting","\n")
      }
      fit <- try(eval(gw.fit(x = X, y = Y, weights = weights, k = k, kstart = fit2$betaIIpars[1], rostart = fit2$betaIIpars[2], betastart = fit2$betascoefs, offset=offset, control = control, method = "L-BFGS-B", hessian = hessian)),silent=TRUE)
      if ('try-error' %in% class(fit)){
        if (control$trace > 0) {
          cat("Crashed 'L-BFGS-B' final fit", "\n")
        }
        fit <- fit2
      }

      fitted <- TRUE

    }
  }
  else{
    if (!any(method == c("nlm", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))){
      stop("method must be in c('nlm', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN')")
    }
    else{
      fit <- try(eval(gw.fit(x = X, y = Y, weights = weights,offset=offset, k = k, kstart = kstart, rostart = rostart, betastart = betastart, control = control, method = method, hessian = hessian)), silent = TRUE)
      if ('try-error' %in% class(fit)){
        warning("Crashed fit", "\n")
      }
      else{
        fitted<-TRUE
      }
    }
  }
  if (fitted){
    if (x)
      fit$X <- X
    if (!y)
      fit$Y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = Terms, data = data, offset = offset, control = control, method = method, contrasts = attr(X, "contrasts"), xlevels = .getXlevels(Terms, mf)))
  }
  else{
    fit<-list()
    fit$aic<- -Inf
    fit$converged=FALSE
    fit$nobs=nrow(X)
  }
  options(warn=warningDefault)
  class(fit) <- "gw"
  fit
}

#' @importFrom stats nlm optim
#' @export
gw.fit <-function (x, y, weights = NULL, k = NULL, kstart = 1, rostart = 2, betastart = NULL, offset = NULL, control = list(), method = "L-BFGS-B", hessian=TRUE, intercept = TRUE){

  control <- do.call("gw.control", control)
  X <- x
  xnames<-dimnames(X)[[2L]]
  Y <- y
  conv <- FALSE
  nobs <- nrow(as.matrix(X))
  ncovars <- ncol(as.matrix(X))
  if (is.null(weights))
    w <- rep(1, nobs)
  else{
    w <- weights
  }
  if (is.null(offset)){
    offset <- rep(0, nobs)
    covoffset <- FALSE
  }
  else{
    covoffset <- TRUE
  }
  if (is.null(kstart))
    kstart <- 1
  if (is.null(rostart))
    rostart <- 2
  if(!is.null(k)){
    kBool <- TRUE   #Boolean about k parameter
    gCorrect <- 1   #To correct df, since there is a parameter less
  }
  else{
    kBool <- FALSE
    gCorrect <- 0
  }

  if (kstart <= 0){
    stop("kstart must be positive")
  }

  if (rostart <= 1){
    stop("rostart must be greater than 1")
  }

  if (!is.logical(hessian)){
    stop("Hessian must be logical")
  }

  if (!any(method == c("nlm", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))){
    stop("method must be in c('nlm', 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'SANN')")
  }

  if (is.null(betastart)){
    betastart <- rep(0, ncovars)
  }

  #Log-likelihood and initial value
  if(!kBool){
    logL<-function(p){
      beta <- p[1:(ncovars)]
      betak <- p[ncovars + 1]
      betaro <- p[ncovars + 2]
      if (method != "L-BFGS-B"){
        k <- exp(betak)
        ro <- 1 + exp(betaro)
      }
      else{
        k <- betak
        ro <- betaro
      }
      mu <- exp(offset + X %*% beta)
      a <- mu * (ro - 1) / k
      gama <- a + k + ro
      -sum(w * (lgamma(a + Y) - lgamma(a) + lgamma(k + Y) - lgamma(k) - lgamma(gama + Y) + lgamma(a + ro) + lgamma(k + ro) - lgamma(ro)))
    }
    if (method != "L-BFGS-B"){
      p0 <- c(betastart, log(kstart), log(rostart - 1))
    }
    else{
      p0 <- c(betastart, kstart, rostart)
    }
  }
  else{
    logL <- function(p){
      beta <- p[1:(ncovars)]
      betaro<-p[ncovars + 1]
      if (method != "L-BFGS-B"){
        ro <- 1 + exp(betaro)
      }
      else{
        ro <- betaro
      }
      mu <- exp(offset + X %*% beta)
      a <- mu * (ro - 1) / k
      gama <- a + k + ro
      #If k=1 there is a simpler expresion
      if (k == 1){
        -sum(w * (lbeta(a + Y, ro + 1) - lbeta(a, ro)))
      }
      else{
        -sum(w * (lgamma(a + Y) - lgamma(a) + lgamma(k + Y) - lgamma(k) - lgamma(gama + Y) + lgamma(a + ro) + lgamma(k + ro) - lgamma(ro)))
      }
    }
    if (method != "L-BFGS-B"){
      p0 <- c(betastart, log(rostart - 1))
    }
    else{
      p0 <- c(betastart, rostart)
    }
  }

  #Optimizing log-likelihood
  if (method == "nlm"){
    fit <- nlm(logL, p = p0, hessian = hessian, iterlim = control$maxit, print.level = control$trace)
    fit$value <- fit$minimum
    fit$par <- fit$estimate
    fit$convergence <- fit$code
    if(fit$convergence<3)
      fit$converged=TRUE
    else
      fit$converged=FALSE
    methodText = "nlm"
  }
  else if (any(method == c("Nelder-Mead", "BFGS", "CG","SANN"))){
    fit <- optim(p0, logL, method = method, hessian = hessian, control = list(maxit = control$maxit, trace = control$trace))
    methodText <- method
    if(fit$convergence==0)
      fit$converged=TRUE
    else
      fit$converged=FALSE
  }
  else if (any(method == c("L-BFGS-B"))){
    #The constraints depends on whether k is known so
    if(!kBool){
      lower <- c(rep(-Inf, ncovars), 0.0000001, 1.0000001)
    }
    else{
      lower <- c(rep(-Inf, ncovars), 1.0000001)
    }
    fit <- optim(p0, logL, method = method, hessian = hessian, lower = lower, control = list(maxit = control$maxit, trace = control$trace))
    methodText <- method
    if(fit$convergence==0)
      fit$converged=TRUE
    else
      fit$converged=FALSE
  }

  if(!kBool){
    if (any(method == c("nlm", "Nelder-Mead", "BFGS", "CG","SANN"))){
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("", c(xnames, "log(k)", "log(ro-1)"))
      betaIIpars <- c(exp(fit$par[ncovars + 1]), 1 + exp(fit$par[ncovars + 2]))
    }
    else{
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("", c(xnames, "k", "ro"))
      betaIIpars <- c(fit$par[ncovars + 1], fit$par[ncovars + 2])
    }
  }
  else{
    if (any(method == c("nlm", "Nelder-Mead", "BFGS", "CG","SANN"))){
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("", c(rep("beta", ncovars), "log(ro-1)"))
      dimnames(coef.table) <- list("", c(xnames, "log(ro-1)"))
      betaIIpars <- c(k, 1 + exp(fit$par[ncovars + 1]))
    }
    else{
      coef.table <- rbind(fit$par, deparse.level = 0)
      dimnames(coef.table) <- list("", c(xnames, "ro"))
      betaIIpars <- c(k, fit$par[ncovars + 1])
    }
  }


  results <- list(
    Y = Y,
    W = w,
    covars = dimnames(X)[[2]],
    nobs = sum(w),
    covoffset = covoffset,
    loglik = -(fit$value + sum(w * lfactorial(Y))),
    aic = 2 * (fit$value + sum(w * lfactorial(Y))) + (length(betastart) + 2 - gCorrect) * 2,
    bic = 2 * (fit$value + sum(w * lfactorial(Y))) + (length(betastart) + 2 - gCorrect) * log(sum(w)),
    df.residual = sum(w) - (length(betastart) + 2 - gCorrect),
    residuals = Y - exp(offset + X %*% fit$par[1:ncovars]),
    coefficients = coef.table,
    betaIIpars = betaIIpars,
    betascoefs = fit$par[1:ncovars],
    fitted.values = exp(offset + X %*% fit$par[1:ncovars]),
    hessian = fit$hessian,
    cov = solve(fit$hessian),
    se = sqrt(diag(solve(fit$hessian))),
    corr = solve(fit$hessian) / (sqrt(diag(solve(fit$hessian))) %o% sqrt(diag(solve(fit$hessian)))),
    code = fit$convergence,
    converged=fit$converged,
    method = methodText,
    k = k,
    kBool = kBool
  )

  class(results) <- "gw"
  return(results)
}

#' @importFrom stats coef naprint
#' @export
print.gw<-function (x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients")
    if (is.character(co <- x$contrasts)) {
      cat("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
    }
    cat(":\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")

  cat("\nDegrees of Freedom:", x$df.residual, "Residual\n")
  if (nzchar(mess <- naprint(x$na.action))) {
    cat("  (", mess, ")\n", sep = "")
  }
  cat("AIC:", format(signif(x$aic, digits)))
  cat("\n")
  if(!x$converged){
    cat("Error of convergence")
  }
  invisible(x)
}

#' @importFrom stats pnorm
#' @export
summary.gw <- function (object, ...){
  df.r <- object$df.residual
  if (is.null(object$covars)) object$covars <- "(Intercept)"
  coef.p <- object$betascoefs
  s.err <- object$se[1:length(object$betascoefs)]
  tvalue <- (object$betascoefs) / (object$se[1:length(object$betascoefs)])
  pvalue <- 2 * pnorm(abs((object$betascoefs) / (object$se[1:length(object$betascoefs)])), lower.tail = FALSE)
  dn <- c("Estimate", "Std. Error")
  coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
  dimnames(coef.table) <- list(object$covars, c(dn, "z value", "Pr(>|z|)"))
  fitted <- cbind(object$loglik, object$aic, object$bic, object$df)
  dimnames(fitted) <- list("", c("log-likelihood", "AIC", "BIC", "df"))
  if (!object$kBool){
    if (any(object$method == c("nlm", "Nelder-Mead", "BFGS", "CG","SANN"))){
      coefk <- object$coefficients[length(object$betascoefs) + 1]
      coefro <- object$coefficients[length(object$betascoefs) + 2]
      k <- object$betaIIpars[1]
      ro <- object$betaIIpars[2]
      #Std. error aproximated by delta method
      SE.k <- object$se[length(object$betascoefs) + 1] * exp(coefk)
      SE.ro <- object$se[length(object$betascoefs) + 2] * exp(coefro)
      betaII <- cbind(c("k", "ro"), format(c(k, ro)), format(c(SE.k, SE.ro)))
      dimnames(betaII) <- list(c("",""), c("par", "Estimate", "Std. Error (Delta method)"))
    }
    else{
      k <- object$betaIIpars[1]
      ro <- object$betaIIpars[2]
      SE.k <- object$se[length(object$betascoefs) + 1]
      SE.ro <- object$se[length(object$betascoefs) + 2]
      betaII <- cbind(c("k", "ro"), format(c(k, ro)), format(c(SE.k, SE.ro)))
      dimnames(betaII) <- list(c("", ""), c("par", "Estimate", "Std. Error"))
    }
  }
  else{
    if (any(object$method == c("nlm", "Nelder-Mead", "BFGS", "CG","SANN"))){
      coefro <- object$coefficients[length(object$betascoefs) + 1]
      k <- 1
      ro <- object$betaIIpars[2]
      #Std. Error aproximated by delta method
      SE.k <- NA
      SE.ro <- object$se[length(object$betascoefs) + 1] * exp(coefro)
      betaII <- cbind(c("k", "ro"), format(c(k, ro)), format(c(SE.k, SE.ro)))
      dimnames(betaII) <- list(c("",""), c("par", "Estimate", "Std. Error"))
    }
    else{
      k <- 1
      ro <- object$betaIIpars[2]
      SE.k <- NA
      SE.ro <- object$se[length(object$betascoefs) + 1]
      betaII <- cbind(c("k", "ro"), format(c(k, ro)), format(c(SE.k, SE.ro)))
      dimnames(betaII) <- list(c("",""), c("par", "Estimate", "Std. Error"))
    }
  }

  keep <- match(c("call", "terms", "deviance", "aic", "contrasts", "df.residual", "na.action"), names(object), 0L)
  ans <- c(object[keep], list(coefficients = coef.table, fitted = fitted, betaII = betaII, method = object$method, convergence = object$code),converged=object$converged)
  class(ans) <- "summary.gw"
  return(ans)
}

#' @importFrom stats coef naprint
#' @export
print.summary.gw <- function (x, digits = max(3, getOption("digits") - 3), ...){
  if (!inherits(x, "summary.gw"))
  stop("'x' must inherit from class %s", dQuote("summary.table"), domain = NA)

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))){
    cat("Coefficients")
    if (is.character(co <- x$contrasts)){
      cat("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
    }
    cat(":\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")
  if (length(x$fitted)){
    cat("\n")
    cat("Fit:\n")
    print.default(format(x$fitted, digits = digits), print.gap = 2, quote = FALSE)
  }
  else cat("No fits\n\n")
  if (length(x$betaII)){
    cat("\n")
    cat("betaII:\n")
    print.default(format(x$betaII, digits = digits), print.gap = 2, quote = FALSE)
  }
  else cat("No betaII\n\n")
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", x$df.residual, "Residual\n")
  cat("\nCode of convergence:", x$convergence, if(!x$converged) "(Error of convergence)" else "", "\n")
  cat("\nMethod:", x$method, "\n")
  if (nzchar(mess <- naprint(x$na.action))){
    cat("  (", mess, ")\n", sep = "")
  }
  invisible(x)
}

gw.control <- function (maxit = 10000, epsilon = 1e-08, trace = FALSE) {
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  list(epsilon = epsilon, maxit = maxit, trace = trace)
}

#' @importFrom stats model.frame
#' @export
model.matrix.gw<-function (object, ...) {
  if (n_match <- match("x", names(object), 0L)) object[[n_match]]
  else {
    data <- model.frame(object, xlev = object$xlevels, ...)
    NextMethod("model.matrix", data = object$data, contrasts.arg = object$contrasts)
  }
}

model.frame.gw <- function (formula, ...) {
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if (length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- quote(stats::model.frame)
    fcall$xlev <- formula$xlevels
    fcall$formula <- terms(formula)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env))
      env <- parent.frame()
    eval(fcall, env, parent.frame())
  }
  else formula$model
}

formula.gw<-function (x, ...) {
  form <- x$formula
  if (!is.null(form)) {
    form <- formula(x$terms)
    environment(form) <- environment(x$formula)
    form
  }
  else formula(x$terms)
}
