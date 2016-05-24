globalVariables(c("post", "intercept", "penalty", "control", "error", "n", "select.Z" , "select.X", "aes", "element_blank"))

#' rlasso: Function for Lasso estimation under homoscedastic and heterosceadstic non-Gaussian
#' disturbances
#'
#' The function estimates the coefficients of a Lasso regression with
#' data-driven penalty under homoscedasticity and heteroscedasticity with non-Gaussian noise and X-dependent or X-independent design. The
#' method of the data-driven penalty can be chosen. The object which is
#' returned is of the S3 class \code{rlasso}.
#'
#' The function estimates the coefficients of a Lasso regression with
#' data-driven penalty under homoscedasticity / heteroscedasticity and non-Gaussian noise. The options \code{homoscedastic} is a logical with \code{FALSE} by default.
#' Moreover, for the calculation of the penalty parameter it can be chosen, if the penalization parameter depends on the  design matrix (\code{X.dependent.lambda=TRUE}) or \code{independent} (default, \code{X.dependent.lambda=FALSE}).
#' The default value of the constant \code{c} is \code{1.1} in the post-Lasso case and \code{0.5} in the Lasso case. 
#'  A \emph{special} option is to set \code{homoscedastic} to \code{none} and to supply a values \code{lambda.start}. Then this value is used as penalty parameter with independent design and heteroscedastic errors to weight the regressors.
#' For details of the
#' implementation of the Algorithm for estimation of the data-driven penalty,
#' in particular the regressor-independent loadings, we refer to Appendix A in
#' Belloni et al. (2012). When the option "none" is chosen for \code{homoscedastic} (together with
#' \code{lambda.start}), lambda is set to \code{lambda.start} and the
#' regressor-independent loadings und heteroscedasticity are used. The options "X-dependent" and
#' "X-independent" under homoscedasticity are described in Belloni et al. (2013). 
# \code{lambda.start} can be component-specific. When used with one of the
# other option, the values are used as starting values.
#'
#' The option \code{post=TRUE} conducts post-lasso estimation, i.e. a refit of
#' the model with the selected variables.
#'
#' @aliases rlasso
#' @param formula an object of class "formula" (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted in the form
#' \code{y~x}
#' @param data an optional data frame, list or environment (or object coercible
#' by as.data.frame to a data frame) containing the variables in the model. If
#' not found in data, the variables are taken from environment(formula),
#' typically the environment from which \code{rlasso} is called.
#' @param post logical. If \code{TRUE}, post-Lasso estimation is conducted.
#' @param intercept logical. If \code{TRUE}, intercept is included which is not
#' penalized.
#' @param penalty list with options for the calculation of the penalty. 
#' \itemize{
#' \item{\code{c} and \code{gamma}}{ constants for the penalty with default \code{c=1.1} and \code{gamma=0.1}}
#' \item{\code{homoscedastic}}{ logical, if homoscedastic errors are considered (default \code{FALSE}). Option \code{none} is described below.}
#' \item{\code{X.dependent.lambda}}{ logical,  \code{TRUE}, if the penalization parameter depends on the the design of the matrix \code{x}. \code{FALSE}, if independent of the design matrix  (default).}
#' \item{\code{numSim}}{ number of simulations for the dependent methods, default=5000}
#' \item{\code{lambda.start}}{ initial penalization value, compulsory for method "none"}
#' }
#' @param control list with control values.
#' \code{numIter} number of iterations for the algorithm for
#' the estimation of the variance and data-driven penalty, ie. loadings,
#' \code{tol} tolerance for improvement of the estimated variances.
#'\code{threshold} is applied to the final estimated lasso
#' coefficients. Absolute values below the threshold are set to zero.
#' @param ... further arguments (only for consistent defintion of methods)
#' @return \code{rlasso} returns an object of class \code{rlasso}. An object of
#' class "rlasso" is a list containing at least the following components:
#' \item{coefficients}{parameter estimates (named vector of coefficients without intercept)}
#' \item{intercept.value}{value of the intercept}
#' \item{index}{index of selected
#' variables (logical vector)} \item{lambda}{data-driven penalty term for each
#' variable, product of lambda0 (the penalization parameter) and the loadings}
#' \item{lambda0}{penalty term} \item{loadings}{loading for each regressor}
#' \item{residuals}{residuals, response minus fitted values} \item{sigma}{root of the variance of
#' the residuals} \item{iter}{number of iterations} \item{call}{function call}
#' \item{options}{options}
#' @references A. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012).
#' Sparse models and methods for optimal instruments with an application to
#' eminent domain. \emph{Econometrica} 80 (6), 2369-2429.
#' @references A. Belloni, V. Chernozhukov and C. Hansen (2013). Inference for
#' high-dimensional sparse econometric models. In Advances in Economics and
#' Econometrics: 10th World Congress, Vol. 3: Econometrics, Cambirdge
#' University Press: Cambridge, 245-295.
#' @keywords Lasso data-driven penalty non-Gaussian heteroscedasticity
#' @export
#' @rdname rlasso
#' @examples
#' set.seed(1)
#' n = 100 #sample size
#' p = 100 # number of variables
#' s = 3 # nubmer of non-zero variables
#' X = matrix(rnorm(n*p), ncol=p)
#' beta = c(rep(3,s), rep(0,p-s))
#' y = 1 + X%*%beta + rnorm(n)
#' lasso.reg = rlasso(y~X, post=TRUE, intercept=TRUE)
#' lasso.reg = rlasso(y~X,post=FALSE)  # use Lasso, not-Post-Lasso
#' print(lasso.reg, all=FALSE)
#' summary(lasso.reg, all=FALSE) # use this option to summarize results
#' yhat.lasso = predict(lasso.reg)   #in-sample prediction
#' Xnew = matrix(rnorm(n*p), ncol=p)  # new X
#' Ynew =  Xnew%*%beta + rnorm(n)  #new Y
#' yhat.lasso.new = predict(lasso.reg, newdata=Xnew)  #out-of-sample prediction
 rlasso <- function(formula, data, post = TRUE, intercept = TRUE, 
                           penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE, lambda.start = NULL, c = 1.1, gamma = .1/log(n)),
                          control = list(numIter = 15, tol = 10^-5, threshold = NULL), ...) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  attr(mt, "intercept") <- 0
  y <- model.response(mf, "numeric")
  n <- length(y)
  x <- model.matrix(mt, mf)
  est <- rlasso.fit(x, y, post = post, intercept = intercept, penalty=penalty,
               control = control)
  est$call <- cl
  return(est)
}

#' @rdname rlasso
#' @export
#' @param y dependent variable (vector, matrix or object can be coerced to matrix)
#' @param x regressors (vector, matrix or object can be coerced to matrix)
rlasso.fit <- function(x, y, post = TRUE, intercept = TRUE,
                           penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE, lambda.start = NULL, c = 1.1, gamma = 0.1),
                           control = list(numIter = 15, tol = 10^-5, threshold = NULL),...) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  if (is.null(colnames(x)))
    colnames(x) <- paste("V", 1:p, sep = "")
  ind.names <- 1:p
  
  # set options to default values if missing
  if (!exists("homoscedastic", where = penalty))  penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty))  penalty$X.dependent.lambda = "FALSE"
  if (!exists("gamma", where = penalty))  penalty$gamma = 0.1/log(n)
  
  if (penalty$homoscedastic=="none" & !exists("lambda.start", where=penalty)) error("lambda.start must be provided!")
  # checking input numIter, tol
  if (!exists("numIter", where = control)) {
    control$numIter = 15
  }
  
  if (!exists("tol", where = control)) {
    control$tol = 10^-5
  }
  
  if (post==FALSE & (!exists("c", where = penalty) | is.na(match("penalty", names(as.list(match.call)))))) {
    penalty$c = 0.5
  }
  
  # Intercept handling and scaling
  if (intercept) {
    meanx <- colMeans(x)
    x <- scale(x, meanx, FALSE)
    mu <- mean(y)
    y <- y - mu
  } else {
    meanx <- rep(0, p)
    mu <- 0
  }
  
  normx <- sqrt(apply(x, 2, var))
  ind <- rep(FALSE, p) #
  
  # variables with low variation are taken out, because normalization is not reliable
  # eps <- 10^-9  # precision for scaling
  #ind <- which(normx < eps)
  #if (length(ind) != 0) {
  #  x <- x[, -ind]
  #  normx <- normx[-ind]
  #  ind.names <- ind.names[-ind]
  #  p <- dim(x)[2]
  #  if (!is.null(penalty$lambda.start)) {
  #    penalty$lambda.start <- penalty$lambda.start[-ind]
  #  }
  #}
  
  #
  
  XX <- crossprod(x)
  Xy <- crossprod(x, y)
  
  startingval <- init_values(x,y)$residuals
  pen <- lambdaCalculation(penalty = penalty, y = startingval, x = x)
  lambda <- pen$lambda
  Ups0 <- Ups1 <- pen$Ups0
  lambda0 <- pen$lambda0
  
  mm <- 1
  s0 <- sqrt(var(y))
  while (mm <= control$numIter) {
    # calculation parameters
    coefTemp <- LassoShooting.fit(x, y, lambda, XX = XX, Xy = Xy)$coefficients
    coefTemp[is.na(coefTemp)] <- 0
    ind1 <- (abs(coefTemp) > 0)
    x1 <- as.matrix(x[, ind1, drop = FALSE])
    if (dim(x1)[2] == 0) {
      if (intercept) {
        intercept.value <- mean(y + mu)
      } else {
        intercept.value <- mean(y)
      }
      est <- list(coefficients = rep(0, p), intercept.value=intercept.value, index = rep(FALSE, p),
                  lambda = lambda, lambda0 = lambda0, loadings = Ups0, residuals = y -
                    mean(y), sigma = var(y), iter = mm, call = match.call(),
                  options = list(post = post, intercept = intercept, ind.scale=ind, 
                                 control = control, mu = mu, meanx = meanx))
      class(est) <- "rlasso"
      return(est)
    }
    
    # refinement variance estimation
    if (post) {
      reg <- lm(y ~ -1 + x1)
      coefT <- coef(reg)
      coefT[is.na(coefT)] <- 0
      e1 <- y - x1 %*% coefT
      coefTemp[ind1] <- coefT
    }
    if (!post) {
      e1 <- y - x1 %*% coefTemp[ind1]
    }
    s1 <- sqrt(var(e1))
    
    # homoscedatic and X-independent
    if (penalty$homoscedastic == TRUE && penalty$X.dependent.lambda == FALSE) {
      Ups1 <- s1*normx
      lambda <- rep(pen$lambda0 * s1, p)
    }
    # homoscedatic and X-dependent
    if (penalty$homoscedastic == TRUE && penalty$X.dependent.lambda == TRUE) {
      Ups1 <- s1*normx
      lambda <- rep(pen$lambda0 * s1, p)
    }
    # heteroscedastic and X-independent
    if (penalty$homoscedastic == FALSE && penalty$X.dependent.lambda == FALSE) {
      Ups1 <- 1/sqrt(n) * sqrt(t(t(e1^2) %*% (x^2)))
      lambda <- pen$lambda0 * Ups1
    }
    
    # heteroscedastic and X-dependent
    if (penalty$homoscedastic == FALSE && penalty$X.dependent.lambda == TRUE) {
      lc <- lambdaCalculation(penalty, y=e1, x=x)
      Ups1 <- lc$Ups0
      lambda <- lc$lambda
    }
    
    
    
    # none
    if (penalty$homoscedastic == "none") {
      if (is.null(penalty$lambda.start)) stop("Argument lambda.start required!")
      Ups1 <- 1/sqrt(n) * sqrt(t(t(e1^2) %*% (x^2)))
      lambda <- pen$lambda0 * Ups1
    }
    
    mm <- mm + 1
    if (abs(s0 - s1) < control$tol) {
      break
    }
    s0 <- s1
  }
  
  if (dim(x1)[2] == 0) {
    coefTemp = NULL
    ind1 <- rep(0, p)
  }
  coefTemp <- as.vector(coefTemp)
  coefTemp[abs(coefTemp) < control$threshold] <- 0
  ind1 <- as.vector(ind1)
  coefTemp <- as.vector(as.vector(coefTemp))
  names(coefTemp) <- names(ind1) <- colnames(x)
  if (intercept) {
    if (is.null(mu)) mu <-0
    if (is.null(meanx))  meanx <-  rep(0, length(coefTemp))  #<- 0
    if (sum(ind)==0) {
      intercept.value <- mu - sum(meanx*coefTemp)
    } else {
      intercept.value <- mu - sum(meanx*coefTemp) #sum(meanx[-ind]*coefTemp)
    }
  } else {
    intercept.value <- NA
  }
  est <- list(coefficients = coefTemp, intercept.value=intercept.value, index = ind1, lambda = lambda,
              lambda0 = lambda0, loadings = Ups1, residuals = as.vector(e1), sigma = s1,
              iter = mm, call = match.call(), options = list(post = post, intercept = intercept,
                                                             control = control, penalty = penalty, ind.scale=ind,
                                                             mu = mu, meanx = meanx))
  class(est) <- "rlasso"
  return(est)
}


################ function lambdaCalculation

#' Function for Calculation of the penalty parameter
#'
#' This function implements different methods for calculation of the penalization parameter \eqn{\lambda}. Further details can be found under \link{rlasso}.
#'
#' @param penalty list with options for the calculation of the penalty. 
#' \itemize{
#' \item{\code{c} and \code{gamma}}{ constants for the penalty with default \code{c=1.1} and \code{gamma=0.1}}
#' \item{\code{homoscedastic}}{ logical, if homoscedastic errors are considered (default \code{FALSE}). Option \code{none} is described below.}
#' \item{\code{X.dependent.lambda}}{ if \code{independent} or \code{dependent} design matrix \code{X} is assumed for calculation of the parameter \eqn{\lambda}}
#' \item{\code{numSim}}{ number of simulations for the X-dependent methods}
#' \item{\code{lambda.start}}{ initial penalization value, compulsory for method "none"}
#' }
#' @param x matrix of regressor variables
#' @param y residual which is used for calculation of the variance or the data-dependent loadings
#' @return The functions returns a list with the penalty \code{lambda} which is the product of \code{lambda0} and \code{Ups0}. \code{Ups0}
#' denotes either the variance (\code{independent} case) or the data-dependent loadings for the regressors. \code{method} gives the selected method for the calculation.
#' @export


lambdaCalculation <- function(penalty = list(homoscedastic = FALSE, X.dependent.lambda = FALSE, lambda.start = NULL, c = 1.1, gamma = 0.1),
                              y = NULL, x = NULL) {
  checkmate::checkChoice(penalty$X.dependent.lambda, c(TRUE, FALSE, NULL))
  checkmate::checkChoice(penalty$homoscedastic, c(TRUE, FALSE, "none"))
  if (!exists("homoscedastic", where = penalty))  penalty$homoscedastic = "FALSE"
  if (!exists("X.dependent.lambda", where = penalty))  penalty$X.dependent.lambda = "FALSE"
  if (!exists("c", where = penalty) & penalty$homoscedastic!="none") {
    penalty$c = 1.1
  }
  if (!exists("gamma", where = penalty) & penalty$homoscedastic!="none") {
    penalty$gamma = 0.1
  }


  # homoscedastic and X-independent
  if (penalty$homoscedastic==TRUE && penalty$X.dependent.lambda == FALSE) {
    p <- dim(x)[2]
    n <- dim(x)[1]
    lambda0 <- 2 * penalty$c * sqrt(n) * qnorm(1 - penalty$gamma/(2 *
                                                                    p))
    Ups0 <- sqrt(var(y))
    lambda <- rep(lambda0 * Ups0, p)
  }

  # homoscedastic and X-dependent
  if (penalty$homoscedastic==TRUE && penalty$X.dependent.lambda == TRUE) {
    if (!exists("numSim", where = penalty)) {
      penalty$numSim = 5000
    }
    p <- dim(x)[2]
    n <- dim(x)[1]
    R <- penalty$numSim
    sim <- vector("numeric", length = R)
    for (l in 1:R) {
      g <- matrix(rep(rnorm(n), each = p), ncol = p, byrow = TRUE)
      sim[l] <- n * max(2 * colMeans(x * g))
    }
    lambda0 <- penalty$c * quantile(sim, probs = 1 - penalty$gamma)
    Ups0 <- sqrt(var(y))
    lambda <- rep(lambda0 * Ups0, p)
  }

  # heteroscedastic and X-independent (was "standard")
  if (penalty$homoscedastic==FALSE && penalty$X.dependent.lambda == FALSE) {
    p <- dim(x)[2]
    n <- dim(x)[1]
    lambda0 <- 2*penalty$c*sqrt(n)*sqrt(2*log(2*p*log(n)/penalty$gamma))
    #lambda0 <- 2 * penalty$c * sqrt(n) * qnorm(1 - penalty$gamma/(2 *
    #                                                                p * 1))  # 1=num endogenous variables
    Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    lambda <- lambda0 * Ups0
  }
  
  # heteroscedastic and X-dependent
  if (penalty$homoscedastic==FALSE && penalty$X.dependent.lambda == TRUE) {
    if (!exists("numSim", where = penalty)) {
      penalty$numSim = 5000
    }
    p <- dim(x)[2]
    n <- dim(x)[1]
    R <- penalty$numSim
    sim <- vector("numeric", length = R)
    lasso.x.y <- rlasso(y ~ x)
    eh <- lasso.x.y$residuals
    ehat <- matrix(rep(eh, each = p), ncol = p, byrow = TRUE) # might be improved by initial estimator or passed through
    for (l in 1:R) {
      g <- matrix(rep(rnorm(n), each = p), ncol = p, byrow = TRUE)
      sim[l] <- n * max(2 * colMeans(x * ehat* g)) # sqrt(n) or n??
    }
    lambda0 <- penalty$c * quantile(sim, probs = 1 - penalty$gamma)
    Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    lambda <- lambda0 * Ups0
  }
  
  
  if (!is.null(penalty$lambda.start)) {
    p <- dim(x)[2]
    if (length(penalty$lambda.start) == 1) {
      lambda.start <- rep(penalty$lambda.start, p)
    }
    lambda <- as.matrix(penalty$lambda.start)
  }

  if (penalty$homoscedastic == "none") {
    if (is.null(penalty$lambda.start) | !exists("lambda.start", where = penalty))
      error("For method \"none\" lambda.start must be provided")
    n <- dim(x)[1]
    lambda0 <- penalty$lambda.start
    Ups0 <- 1/sqrt(n) * sqrt(t(t(y^2) %*% (x^2)))
    lambda <- lambda0 * Ups0
  }

  return(list(lambda0 = lambda0, lambda = lambda, Ups0 = Ups0, method = penalty))
}





################# Methods for Lasso

#' Methods for S3 object \code{rlasso}
#'
#' Objects of class \code{rlasso} are constructed by \code{rlasso.formula} or \code{rlasso.fit}.
#' \code{print.rlasso} prints and displays some information about fitted \code{rlasso} objects.
#' \code{summary.rlasso} summarizes information of a fitted \code{rlasso} object.
#' \code{predict.rlasso} predicts values based on a \code{rlasso} object.
#' \code{model.matrix.rlasso} constructs the model matrix of a \code{rlasso} object.
#'
#' @param object an object of class \code{rlasso}
#' @param x an object of class \code{rlasso}
#' @param all logical, indicates if coefficients of all variables (TRUE) should be displayed or only the non-zero ones (FALSE)
#' @param digits significant digits in printout
#' @param newdata new data set for prediction. An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are returned.
#' @param ... arguments passed to the print function and other methods
#' @keywords methods rlasso
#' @rdname methods.rlasso
#' @aliases methods.rlasso print.rlasso predict.rlasso model.matrix.rlasso
#' @export

print.rlasso <- function(x, all=TRUE ,digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    if (all) {
      cat("Coefficients:\n")
      print.default(format(coef(x), digits = digits), print.gap = 2L,
                    quote = FALSE)
    } else {
      print.default(format(coef(x)[x$index], digits = digits), print.gap = 2L,
                    quote = FALSE)
    }
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

#' @rdname methods.rlasso
#' @export

summary.rlasso <- function(object, all=TRUE, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nPost-Lasso Estimation: ",  paste(deparse(object$options$post), sep = "\n", collapse = "\n"), "\n", sep = " ")
  coefs <- object$coefficients
  p <- length(coefs)
  num.selected <- sum(abs(object$coefficients)>0)
  cat("\nTotal number of variables:", p)
  cat("\nNumber of selected variables:", num.selected, "\n", sep=" ")
  resid <- object$residuals
  cat("\nResiduals: \n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure(apply(t(resid), 1L, quantile), dimnames = list(nam, dimnames(resid)[[2L]]))
  print(drop(t(rq)), digits = digits)
  cat("\n")
  if (all) {
    coefm <- matrix(NA, length(coefs), 1)
    coefm[,1] <- coefs
    colnames(coefm) <- "Estimate"
    rownames(coefm) <- names(coefs)
    printCoefmat(coefm, digits = digits, na.print = "NA")
  } else {
    coefs <- coefs[abs(coefs)>0]
    coefm <- matrix(NA, length(coefs), 1)
    coefm[,1] <- coefs
    colnames(coefm) <- "Estimate"
    rownames(coefm) <- names(coefs)
    printCoefmat(coefm, digits = digits, na.print = "NA")
  }
  cat("\nResidual standard error:", format(signif(object$sigma, digits)))
  cat("\n")
  invisible(object)
}

#' @rdname methods.rlasso
#' @export

model.matrix.rlasso <- function(object, ...) {
  if (is.call(object$call[[2]])) {
    if(is.null(object$call$data)){
      X <- model.frame(object$call[[2]])
      mt <- attr(X, "terms")
      attr(mt, "intercept") <- 0
      mm <- model.matrix(mt, X)
    } else {
      dataev <- eval(object$call$data)
      mm <- as.matrix(dataev[,names(object$coefficients)])
    }
  } else {
    mm <- eval(object$call[[2]])
  }
  return(mm)
}


#' @rdname methods.rlasso
#' @export

predict.rlasso <- function (object, newdata = NULL, ...){
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)

    if(sum(object$options$ind.scale)!=0) {
      X <- X[,-object$options$ind.scale]
    }
  }
  else {
    varcoef <- names(object$coefficients)
    if(all(is.element(varcoef,colnames(newdata)))){
      X <- as.matrix(newdata[,varcoef])
    } else {
      X <- as.matrix(newdata)

      if(sum(object$options$ind.scale)!=0) {
        X <- X[,-object$options$ind.scale]
      }
      stopifnot(ncol(X)==length(object$coefficients))
    }
  }
  n <- dim(X)[1] #length(object$residuals)
  beta <- object$coefficients

  if (object$options[["intercept"]]) {
    yhat <- X %*% beta + object$intercept.value
    if (dim(X)[2]==0) yhat <- rep(object$intercept.value, n)
  }
  if (!object$options[["intercept"]]) {
    yhat <- X %*% beta
    if (dim(X)[2]==0) yhat <- rep(0, n)
  }
  return(yhat)
}












