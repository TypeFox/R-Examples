#' @title Auxiliary for Controlling ZIM Fitting
#' @description Auxiliary function for \code{\link{zim}} fitting. Typically only used internally by 
#' \code{\link{zim.fit}}, but may be used to construct a control argument for either function.
#' @param dist count model family.
#' @param method algorithm for parameter estimation.
#' @param type type of matrix inverse.
#' @param robust logical; if TRUE, robust standard errors will be calculated.
#' @param trace logical; if TRUE, display iteration history.
#' @param start initial parameter values.
#' @param minit minimum number of iterations.
#' @param maxit maximum number of iterations.
#' @param epsilon positive convergence tolerance.
#' @seealso \code{\link{zim}}, \code{\link{zim.fit}}
#' @keywords regression
#' @export zim.control
zim.control <- function(dist = c("zip", "zinb"), method = c("EM-NR", "EM-FS"), 
  type = c("solve", "ginv"), robust = FALSE, trace = FALSE, start = NULL, 
  minit = 10, maxit = 10000, epsilon = 1e-8) {
  dist <- match.arg(dist)
  method <- match.arg(method)
  type <- match.arg(type)
  inv <- switch(type, "solve" = solve, "ginv" = ginv)
  list(dist = dist, method = method, type = type, inv = inv, robust = robust, 
    trace = trace, start = start, minit = minit, maxit = maxit, epsilon = epsilon)
}

#' @title Fitter Function for Zero-Inflated Models
#' @description \code{\link{zim.fit}} is the basic computing engine called by \code{\link{zim}} used to fit 
#' zero-inflated models. This should usually \emph{not} be used directly unless by experienced users. 
#' @param y response variable.
#' @param X design matrix for log-linear part.
#' @param Z design matrix for logistic part.
#' @param weights an optional vector of 'prior weights' to be used in the fitting process.
#' @param offset offset variable
#' @param control control arguments from \code{\link{zim.control}}.
#' @param ... additional argumetns.
#' @seealso \code{\link{zim}}, \code{\link{zim.control}}
#' @keywords regression
#' @export zim.fit
zim.fit <- function(y, X, Z, weights = rep(1, nobs), offset = rep(0, nobs), 
  control = zim.control(...), ...) { 
  nobs <- length(y)
  pX <- NCOL(X)
  pZ <- NCOL(Z)
  y0 <- (y == 0)
  if(is.null(weights)) weights <- 1
  if(is.null(offset)) offset <- 0
  if(length(weights) == 1) weights <- rep(weights, nobs) 
  if(length(offset) == 1) offset <- rep(offset, nobs)
  loglik <- function(para) {
    zip.loglik <- function(para) {
      lambda <- exp(X %*% para[1:pX] + offset)
      omega <- plogis(Z %*% para[pX + 1:pZ])
      sum(weights * dzip(y, lambda, omega, log = TRUE))
    }
    zinb.loglik <- function(para) {
      k <- exp(para[1])
      lambda <- exp(X %*% para[1 + 1:pX] + offset)
      omega <- plogis(Z %*% para[1 + pX + 1:pZ]) 
      sum(weights * dzinb(y, k, lambda, omega, log = TRUE))
    }
    switch(control$dist, "zip" = zip.loglik(para), "zinb" = zinb.loglik(para))
  }
  em <- function(para) {
    zip.em <- function(para) { 
      lambda <- exp(X %*% para[1:pX] + offset)
      omega <- plogis(Z %*% para[pX + 1:pZ]) 
      p0 <- dzip(0, lambda, omega)
      u <- y0 * omega / p0
      para[1:pX] <- suppressWarnings(glm.fit(X, y, weights = weights * (1 - u), 
        offset = offset, family = poisson(), start = para[1:pX]))$coef
      para[pX + 1:pZ] <- suppressWarnings(glm.fit(Z, u, weights = weights, 
        family = binomial(), start = para[pX + 1:pZ]))$coef
      para
    }  
    zinb.em <- function(para) {
      k <- exp(para[1])
      lambda <- exp(X %*% para[1 + 1:pX] + offset)
      omega <- plogis(Z %*% para[1 + pX + 1:pZ])
      p0 <- dzinb(0, k, lambda, omega)
      u <- y0 * omega / p0
      nb.fit <- suppressWarnings(glm.nb(y ~ 0 + X + offset(offset),
        weights = weights * (1 - u), init.theta = k, start = para[1 + 1:pX]))
      para[1] <- log(nb.fit$theta)
      para[1 + 1:pX] <- nb.fit$coef
      para[1 + pX + 1:pZ] <- suppressWarnings(glm.fit(Z, u, weights = weights, 
        family = binomial(), start = para[1 + pX + 1:pZ]))$coef
      para
    }  
    switch(control$dist, "zip" = zip.em(para), "zinb" = zinb.em(para))
  }
  deriv <- function(para) {
    zip.deriv <- function(para) { 
      lambda <- exp(X %*% para[1:pX] + offset)
      omega <- plogis(Z %*% para[pX + 1:pZ]) 
      p0 <- dzip(0, lambda, omega)    
      v1 <- y - lambda * (1 - omega * y0 / p0)
      v2 <- omega * (y0 / p0 - 1)
      gradient <- rep(NA, pX  + pZ)     
      gradient[1:pX] <- t(X) %*% (weights * v1)
      gradient[pX + 1:pZ] <- t(Z) %*% (weights * v2)     
      J <- matrix(NA, pX + pZ, pX + pZ)
      J[1:pX, 1:pX] <- t(X) %*% (weights * as.vector(v1 * v1) * X)
      J[1:pX, pX + 1:pZ] <- t(X) %*% (weights * as.vector(v1 * v2) * Z)
      J[pX + 1:pZ, pX + 1:pZ] <- t(Z) %*% (weights * as.vector(v2 * v2) * Z)
      J[pX + 1:pZ, 1:pX] <- t(J[1:pX, pX + 1:pZ])
      J <- (J + t(J)) / 2
      info <- matrix(NA, pX + pZ, pX + pZ)
      if(control$method == "EM-NR") {
        info[1:pX, 1:pX] <- t(X) %*% (weights * as.vector(lambda * (1 - y0 * 
          omega * (omega + (1 - omega) * (1 + lambda) * exp(-lambda)) / p0^2)) * X)
        info[1:pX, pX + 1:pZ] <- t(X) %*% (weights * as.vector((-y0) * omega * 
          (1 - omega) * lambda * exp(-lambda) / p0^2) * Z)
        info[pX + 1:pZ, pX + 1:pZ] <- t(Z) %*% (weights * as.vector(omega * 
          (1 - omega) * (1 - y0 * exp(-lambda) / p0^2)) * Z)
        info[pX + 1:pZ, 1:pX] <- t(info[1:pX, pX + 1:pZ])     
      } else {
        info[1:pX, 1:pX] <- t(weights * as.vector((1 - omega) * lambda * 
          (exp(-lambda) + omega * (1 - (1 + lambda) * 
          exp(-lambda))) / p0) * X) %*% X
        info[1:pX, pX + 1:pZ] <- t(weights * as.vector((-1) * omega * 
          (1 - omega) * lambda * exp(-lambda) / p0) * X) %*% Z
        info[pX + 1:pZ, pX + 1:pZ] <- t(weights * as.vector(omega^{2} * 
          (1 - omega) * (1 - exp(-lambda)) / p0) * Z) %*% Z 
        info[pX + 1:pZ, 1:pX] <- t(info[1:pX, pX + 1:pZ])     
      }                                    
      info <- (info + t(info)) / 2      
      list(gradient = gradient, J = J, info = info)
    }  
    zinb.deriv <- function(para) {
      k <- exp(para[1])
      lambda <- exp(X %*% para[1 + 1:pX] + offset)
      omega <- plogis(Z %*% para[1 + pX + 1:pZ])
      p <- k / (k + lambda)
      p0 <- dzinb(0, k, lambda, omega)
      v1 <- k * (1 - omega * y0 / p0) * (log(p) + 1 - 
        p) + k * (1 - y0) * (digamma(k + y) - digamma(k)) - p * y
      v2 <- p * y - k * (1 - p) * (1 - omega * y0 / p0)
      v3 <- omega * (y0 / p0 - 1)
      gradient <- rep(NA, 1 + pX + pZ)
      gradient[1] <- sum(weights * v1)
      gradient[1 + 1:pX] <- t(X) %*% (weights * v2)                                                                          
      gradient[1 + pX + 1:pZ] <- t(Z) %*% (weights * v3)
      J <- matrix(NA, 1 + pX + pZ, 1 + pX + pZ)
      J[1, 1] <- sum(weights * as.vector(v1 * v1))
      J[1, 1 + 1:pX] <- t(X) %*% (weights * as.vector(v1 * v2))
      J[1, 1 + pX + 1:pZ] <- t(Z) %*% (weights * as.vector(v1 * v3))
      J[1 + 1:pX, 1 + 1:pX] <- t(X) %*% (weights * as.vector(v2 * v2) * X)
      J[1 + 1:pX, 1 + pX + 1:pZ] <- t(X) %*% (weights * as.vector(v2 * v3) * Z)
      J[1 + pX + 1:pZ, 1 + pX + 1:pZ] <- t(Z) %*% (weights * as.vector(v3 * v3) * Z)
      J[1 + 1:pX, 1] <- t(J[1, 1 + 1:pX])
      J[1 + pX + 1:pZ, 1] <- t(J[1, 1 + pX + 1:pZ])
      J[1 + pX + 1:pZ, 1 + 1:pX] <- t(J[1 + 1:pX, 1 + pX + 1:pZ])
      J <- (J + t(J)) / 2
      info <- matrix(NA, 1 + pX + pZ, 1 + pX + pZ)
      if(control$method == "EM-NR") {
        info[1, 1] <- sum(weights * (k * (1 - p) * (y * p / k - (1 - 
          p) * (1 - omega * y0 / p0)) - k * (1 - y0) * (digamma(k + y) + 
          k * trigamma(k + y) - digamma(k) - k * trigamma(k)) - 
          k * (log(p) + 1 - p) * (1 - omega * y0 / p0 + k * omega * y0 * 
          (log(p) + 1 - p) * (p0 - omega) / p0^2)))
        info[1, 1 + 1:pX] <- t(X) %*% (weights * (k * (1 - p) * 
          (k * omega * y0 * (log(p) + 1 - p) * (p0 - omega) / p0^2 + 
          (1 - p) * (1 - omega * y0 / p0) - y * p / k)))
        info[1, 1 + pX + 1:pZ] <- t(Z) %*% (weights * (k * omega * y0 * 
          (log(p) + 1 - p) * (p0 - omega) / p0^2))
        info[1 + 1:pX, 1 + 1:pX] <- t(X) %*% (weights * as.vector(k * 
          (1 - p) * (p * (1 + y / k - omega * y0 / p0) - 
          k * omega * y0 * (1 - p) * (p0 - omega) / p0^2)) * X)
        info[1 + 1:pX, 1 + pX + 1:pZ] <- t(X) %*% (weights * 
          as.vector((-k) * omega * (1 - omega) * 
          (1 - p) * p^k * y0 / p0^2) * Z)      
        info[1 + pX + 1:pZ, 1 + pX + 1:pZ] <- t(Z) %*% (weights * 
          as.vector(omega * (1 - omega) * (1 - p^k * y0 / p0^2)) * Z)
        info[1 + 1:pX, 1] <- t(info[1, 1 + 1:pX])
        info[1 + pX + 1:pZ, 1] <- t(info[1, 1 + pX + 1:pZ])
        info[1 + pX + 1:pZ, 1 + 1:pX] <- t(info[1 + 1:pX, 1 + pX + 1:pZ])    
      } else {
        f <- function(v) {
          k <- v[1]
          lambda <- v[2]
          omega <- v[3]
          j <- 0:100
          sum(pzinb(j, k, lambda, omega, lower.tail = FALSE) / (k + j)^2)
        }
        c <- apply(cbind(k, lambda, omega), MARGIN = 1, FUN = f)
        info[1, 1] <- sum(weights * (k * (k * c - (1 - omega) * (1 - p)) - 
          omega * (1 - omega) * k^2 * (log(p) + 1 - p)^2 * p^k / p0))
        info[1, 1 + 1:pX] <- t(X) %*% (weights * (omega * (1 - omega) * 
          k^2 * (1 - p) * (log(p) + 1 - p) * p^k / p0))
        info[1, 1 + pX + 1:pZ] <- t(Z) %*% (weights * (omega * (1 - omega) * 
          k * (log(p) + 1 - p) * p^k / p0))
        info[1 + 1:pX, 1 + 1:pX] <- t(X) %*% (weights * as.vector((1 - omega) * 
          k * (1 - p) * (p^k + omega * (1 - p^k - k * (1 - p) * p^k)) / p0) * X)
        info[1 + 1:pX, 1 + pX + 1:pZ] <- t(X) %*% (weights * 
          as.vector(-omega * (1 - omega) * k * (1 - p) * p^k / p0) * Z)      
        info[1 + pX + 1:pZ, 1 + pX + 1:pZ] <- t(Z) %*% (weights * 
          as.vector(omega^2 * (1 - omega) * (1 - p^k) / p0) * Z)
        info[1 + 1:pX, 1] <- t(info[1, 1 + 1:pX])
        info[1 + pX + 1:pZ, 1] <- t(info[1, 1 + pX + 1:pZ])
        info[1 + pX + 1:pZ, 1 + 1:pX] <- t(info[1 + 1:pX, 1 + pX + 1:pZ])
      }
      info <- (info + t(info)) / 2  
      list(gradient = gradient, J = J, info = info)
    }  
    switch(control$dist, "zip" = zip.deriv(para), "zinb" = zinb.deriv(para))    
  }
  iter <- 0
  if(is.null(control$start)) {
    if(control$dist == "zip") {
      para.old <- rep(NA, pX + pZ)
      para.old[1:pX] <- suppressWarnings(glm.fit(X, y, weights = weights,
        offset = offset, family = poisson()))$coef
      para.old[pX + 1:pZ] <- suppressWarnings(glm.fit(Z, y0, weights = weights,
        family = binomial()))$coef     
    } else {
      para.old <- rep(NA, 1 + pX + pZ)
      nb.fit <- suppressWarnings(glm.nb(y ~ 0 + X + offset(offset),  
        weights = weights))
      para.old[1] <- log(nb.fit$theta)
      para.old[1 + 1:pX] <- nb.fit$coef
      para.old[1 + pX + 1:pZ] <- suppressWarnings(glm.fit(Z, y0, 
        weights = weights, family = binomial()))$coef    
    }
  } else {
    para.old <- control$start
  }
  repeat{
    iter <- iter + 1
    if(iter > control$maxit) {
      stop("The maximum number of iterations has been reached!")  
    }
    if(iter <= control$minit) {
      para.new <- em(para.old)
      loglik.new <- loglik(para.new)  
      if(control$trace) {
        cat("iter =", iter, "\t loglik =", loglik.new, "\n")    
        print(para.new)
        cat("\n")
      }
      para.old <- para.new
    } else {
      deriv.old <- deriv(para.old)
      para.new <- as.vector(para.old +
        control$inv(deriv.old$info) %*% deriv.old$gradient)
      loglik.old <- loglik(para.old)
      loglik.new <- loglik(para.new)
      if(control$trace) {
        cat("iter =", iter, "\t loglik =", loglik.new, "\n")    
        print(para.new)
        cat("\n")
      }
      if(abs(loglik.new - loglik.old) / (abs(loglik.new) + control$epsilon) 
        > control$epsilon) {
        para.old <- para.new 
      } else break         
    }  
  }
  deriv.new <- deriv(para.new) 
  dim <- (control$dist == "zinb") + pX + pZ
  aic <- (-2) * loglik.new + 2 * dim
  bic <- (-2) * loglik.new + log(sum(weights)) * dim
  tic <- (-2) * loglik.new + 2 * sum(diag(deriv.new$J %*% control$inv(deriv.new$info)))
  I.inv <- control$inv(deriv.new$info)
  if(control$robust == TRUE) {
    se <- sqrt(diag(I.inv %*% deriv.new$J %*% I.inv)) 
  } else {
    se <- sqrt(diag(I.inv))
  }
  list(iter = iter, gradient = deriv.new$gradient, info = deriv.new$info,
    para = para.new, se = se, loglik = loglik.new, aic = aic, bic = bic, tic = tic)  
}

#' @title Fitting Zero-Inflated Models
#' @description \code{zim} is used to fit zero-inflated models.
#' @param formula an objective of class "\code{\link{formula}}".
#' @param data an optional dataframe, list or environment containing the variables in the model.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s.
#' @param weights an optional vector of 'prior weights' to be used in the fitting process.
#' @param offset this can be used to specify a priori known component to be included in the linear predictor during fitting.
#' @param control control arguments.
#' @param ... additional arguments.
#' @note \code{\link{zim}} is very similar to \code{\link[pscl]{zeroinfl}} from the \code{pscl} package. Both functions can be used to 
#' fit observation-driven models for zero-inflated time series.  
#' @seealso \code{\link{zim.fit}}, \code{\link{zim.control}}
#' @keywords regression
#' @export zim
zim <- function(formula, data, subset, na.action, weights = 1, offset = 0, 
  control = zim.control(...),  ...) {
  call <- match.call()
  if(missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action",
    "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf[[1]] <- as.name("model.frame")
  mf$drop.unused.levels <- TRUE 
  if(any(as.character(formula[[3]]) == "|")) {
    formulaX <- . ~ .
    formulaZ <- . ~ .
    formulaX[[2]] <- formula[[2]]
    formulaZ[[2]] <- formula[[2]]
    formulaX[[3]] <- formula[[3]][[2]]
    formulaZ[[3]] <- formula[[3]][[3]]
    formula[[3]][[1]] <- as.name("+")
    mf$formula <- formula
  } else {
    formulaX <- formula
    formulaZ <- formula
    mf$formula <- formula
  }
  mf <- eval(mf, parent.frame())    
  y <- round(model.response(mf, "numeric"))
  X <- model.matrix(terms(formulaX, data = data), mf)
  Z <- model.matrix(terms(formulaZ, data = data), mf)
  nobs <- length(y)
  pX <- NCOL(X)
  pZ <- NCOL(Z)
  y0 <- (y == 0)
  weights <- as.vector(model.weights(mf))
  offset <- as.vector(model.offset(mf))
  if(is.null(weights)) {
    weights <- rep(1, nobs)
  }
  if(is.null(offset)) {
    offset <- rep(0, nobs)
  }    
  fit <- zim.fit(y, X, Z, weights = weights, offset = offset, control = control)
  fit$call <- call
  fit$control <- control
  fit$na.action <- attr(mf, "na.action")
  fit$y <- y
  fit$X <- X
  fit$Z <- Z                         
  if(control$dist == "zip") {
    lambda <- exp(X %*% fit$para[1:pX] + offset)
    omega <- plogis(Z %*% fit$para[pX + 1:pZ]) 
    p0 <- dzip(0, lambda, omega)
    J <- matrix(NA, 1 + pX + pZ, 1 + pX + pZ)
    J[1, 1] <- sum(weights * (lambda^2 * (2 * (1 - omega) -
      omega * lambda^2 * (1 - omega / p0)))) / 4
    J[1, 1 + 1:pX] <- t(X) %*% (weights * (omega *
      lambda^3 * (1 - omega / p0))) / 2
    J[1, 1 + pX + 1:pZ] <- t(Z) %*% (weights * (omega *
      lambda^2 * (1 - omega / p0))) / 2
    J[1 + 1:pX, 1 + 1:pX] <- t(X) %*% (weights * as.vector(lambda *
      ((1 - omega) - omega * lambda * (1 - omega / p0))) * X)
    J[1 + 1:pX, 1 + pX + 1:pZ] <- (-1) * t(X) %*% (weights *
      as.vector(omega * lambda * (1 - omega / p0)) * Z)
    J[1 + pX + 1:pZ, 1 + pX + 1:pZ] <- t(Z) %*% (weights *
      as.vector(omega^2 * (1 / p0 - 1)) * Z)
    J[1 + 1:pX, 1] <- t(J[1, 1 + 1:pX])
    J[1 + pX + 1:pZ, 1] <- t(J[1, 1 + pX + 1:pZ])
    J[1 + pX + 1:pZ, 1 + 1:pX] <- t(J[1 + 1:pX, 1 + pX + 1:pZ])
    J <- (J + t(J)) / 2
    S <- sum(weights * ((y - lambda)^2 - y - y0 * omega * lambda^2 / p0)) / 2  
    fit$score.test <- S * sqrt(control$inv(J)[1, 1])
    fit$p.value <- pnorm(fit$score.test, lower.tail = FALSE)
    fit$k <- Inf
    fit$lambda <- lambda
    fit$omega <- omega     
  } else {
    fit$k <- exp(fit$para[1])
    fit$lambda <- exp(X %*% fit$para[1 + 1:pX] + offset)
    fit$omega <- plogis(Z %*% fit$para[1 + pX + 1:pZ]) 
  }
  fit$fitted.values <- fit$lambda * (1 - fit$omega)
  fit$residuals <- (y - fit$fitted.values) / sqrt(fit$fitted.values * 
    (1 + fit$omega * fit$lambda + fit$lambda / fit$k))             
  class(fit) <- "zim"
  fit
}

#' @export
print.zim <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  pX <- NCOL(x$X)
  pZ <- NCOL(x$Z)
  z.value <- x$para / x$se 
  z.prob <- pvalue(z.value)
  coef <- data.frame(x$para, x$se, z.value, z.prob)
  coefX <- coef[(x$control$dist == "zinb") + 1:pX, ]
  coefZ <- coef[(x$control$dist == "zinb") + pX + 1:pZ, ]
  colnames(coefX) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)") 
  colnames(coefZ) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")  
  rownames(coefX) <- colnames(x$X)
  rownames(coefZ) <- colnames(x$Z)     
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")  
  cat("Coefficients (log-linear): \n")   
  printCoefmat(coefX, signif.legend = FALSE)
  cat("\n")
  cat("Coefficients (logistic): \n")     
  printCoefmat(coefZ, signif.legend = FALSE)
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n\n")  
  if(x$control$dist == "zip") {                   
    cat("Test for overdispersion (H0: ZIP vs. H1: ZINB) \n")
    cat("score.test:", round(x$score.test, 4), "\n")
    cat("p.value:", format.pval(x$p.value), "\n\n")
  } else {
    cat(paste("(Dispersion parameter for negative binomial taken to be ",
      round(x$k, 4), ")", sep = ""), "\n\n")
  }
  cat("Criteria for assessing goodness of fit \n") 
  cat("loglik:", x$loglik, "\n")
  cat("aic:", x$aic, "\n")
  cat("bic:", x$bic, "\n") 
  cat("tic:", x$tic, "\n") 
  cat("\n")
  cat("Number of", x$control$method, "iterations:", x$iter, "\n")
  cat("Maximum absolute gradient:", max(abs(x$gradient)), "\n")
  cat("\n")
  invisible(x)
}

