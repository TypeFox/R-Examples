
#' rlassologit: Function for logistic Lasso estimation
#'
#' The function estimates the coefficients of a logistic Lasso regression with
#' data-driven penalty. The method of the data-driven penalty can be chosen.
#' The object which is returned is of the S3 class \code{rlassologit}
#'
#' The function estimates the coefficients of a Logistic Lasso regression with
#' data-driven penalty. The
#' option \code{post=TRUE} conducts post-lasso estimation, i.e. a refit of the
#' model with the selected variables.
#' @param formula an object of class 'formula' (or one that can be coerced to
#' that class): a symbolic description of the model to be fitted in the form
#' \code{y~x}
#' @param data an optional data frame, list or environment
#' @param x regressors (matrix)
#' @param y dependent variable (vector or matrix)
#' @param post logical. If \code{TRUE}, post-lasso estimation is conducted.
#' @param intercept logical. If \code{TRUE}, intercept is included which is not
#' penalized.
#' @param penalty list with options for the calculation of the penalty.  \code{c} and \code{gamma} constants for the penalty.
#' @param control list with control values.
#' \code{threshold} is applied to the final estimated lasso
#' coefficients. Absolute values below the threshold are set to zero.
#' @param ... further parameters passed to glmnet
#' @return \code{rlassologit} returns an object of class
#' \code{rlassologit}. An object of class \code{rlassologit} is a list
#' containing at least the following components: \item{coefficients}{parameter
#' estimates (without intercept)} \item{a0}{value of intercept} \item{index}{index of selected variables (logicals)}
#' \item{lambda}{penalty term}
#' \item{residuals}{residuals}
#' \item{sigma}{root of the variance of the residuals}
#' \item{call}{function call}
#' \item{options}{options}
#' @references Belloni, A., Chernozhukov and Y. Wei (2013). Honest confidence regions for logistic regression with a large number of controls. arXiv preprint arXiv:1304.3969.
#' @keywords logistic lasso lasso logistic regression
#' @export
#' @rdname rlassologit
#' @examples
#'\dontrun{
#' library(hdm)
#' ## DGP
#' set.seed(2)
#' n <- 250
#' p <- 100
#' px <- 10
#' X <- matrix(rnorm(n*p), ncol=p)
#' beta <- c(rep(2,px), rep(0,p-px))
#' intercept <- 1
#' P <- exp(intercept + X %*% beta)/(1+exp(intercept + X %*% beta))
#' y <- numeric(length=250)
#' for(i in 1:n){
#'   y[i] <- sample(x=c(1,0), size=1, prob=c(P[i],1-P[i]))
#' }
#' ## fit rlassologit object
#'  rlassologit.reg <- rlassologit(y~X)
#'  ## methods
#' summary(rlassologit.reg, all=F)
#' print(rlassologit.reg)
#' predict(rlassologit.reg, type='response')
#' X3 <- matrix(rnorm(n*p), ncol=p)
#' predict(rlassologit.reg, newdata=X3)
#' }
#' @export
rlassologit <- function(formula, data, post = TRUE, intercept = TRUE, penalty = list(lambda = NULL, 
                                                                                     c = 1.1, gamma = 0.1/log(n)), control = list(threshold = NULL), ...) {
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
  
  est <- rlassologit.fit(x, y, post = post, intercept = intercept, penalty = penalty, 
                         control = control)
  est$call <- cl
  return(est)
}


#' @rdname rlassologit
#' @export

rlassologit.fit <- function(x, y, post = TRUE, intercept = TRUE, penalty = list(lambda = NULL, 
                                                                                c = 1.1, gamma =  0.1/log(n)), control = list(threshold = NULL), ...) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (is.null(colnames(x))) 
    colnames(x) <- paste("V", 1:p, sep = "")
  ind.names <- 1:p
  
  if (!exists("c", where = penalty)) {
    if (post==TRUE) {
      penalty$c = 1.1
    } else {
      penalty$c = 0.5
    }
  }
  
  if (!exists("gamma", where = penalty)) {
    penalty$gamma = 0.1/log(n)
  }
  
  if (is.null(penalty$gamma)) {
    penalty$gamma = 0.1/log(n)
  }
  
  if (!is.null(penalty$lambda)) {
    lambda <- penalty$lambda/(2 * n)
    lambda0 <- lambda * (2 * n)
  } else {
    lambda0 <- penalty$c/2 * sqrt(n) * qnorm(1 - penalty$gamma/(2 * 
                                                                  p))
    lambda <- lambda0/(2 * n)
  }
  
  s0 <- sqrt(var(y))
  # calculation parameters
  xs <- scale(x, center = FALSE, scale = TRUE)
  log.lasso <- glmnet::glmnet(xs, y, family = c("binomial"), alpha = 1, 
                              lambda = lambda[1], standardize = TRUE, intercept = intercept)
  coefTemp <- as.vector(log.lasso$beta)
  coefTemp[is.na(coefTemp)] <- 0
  ind1 <- (abs(coefTemp) > 0)
  x1 <- as.matrix(x[, ind1, drop = FALSE])
  if (dim(x1)[2] == 0) {
    if (intercept == TRUE) {
      a0 <- log(mean(y)/(1 - mean(y)))
      res <- y - mean(y)
    }
    
    if (intercept == FALSE) {
      a0 <- 0  # or NA?
      res <- rep(y - 0.5, n)
      message("Residuals not defined, set to 0.5")
    }
    est <- list(coefficients = rep(0, p), a0 = a0, index = rep(FALSE, 
                                                               p), s0 = s0, lambda0 = lambda0, residuals = res, sigma = sqrt(var(res)), 
                call = match.call(), options = list(post = post, intercept = intercept, 
                                                    control = control))
    class(est) <- c("rlassologit")
    return(est)
  }
  
  # refinement variance estimation
  if (post) {
    if (intercept) {
      # reg <- glm.fit(x1,y, family = binomial(link = 'logit'),
      # intercept=intercept)
      reg <- glm(y ~ x1, family = binomial(link = "logit"))
      coefT <- coef(reg)[-1]
      coefT[is.na(coefT)] <- 0
      e1 <- y - reg$fitted.values
      coefTemp[ind1] <- coefT
    }
    
    if (!intercept) {
      reg <- glm(y ~ -1 + x1, family = binomial(link = "logit"))
      coefT <- coef(reg)
      coefT[is.na(coefT)] <- 0
      e1 <- y - reg$fitted.values
      coefTemp[ind1] <- coefT
    }
  }
  if (!post) {
    e1 <- y - predict(log.lasso, newx = x, type = "response")
  }
  
  
  if (intercept == TRUE) {
    if (post == TRUE) 
      a0 <- coef(reg)[1]
    if (post == FALSE) 
      a0 <- coef(log.lasso)[1]
  }
  
  if (intercept == FALSE) {
    a0 <- 0  # or NA?
  }
  
  coefTemp <- as.vector(coefTemp)
  coefTemp[abs(coefTemp) < control$threshold] <- 0
  ind1 <- as.vector(ind1)
  names(coefTemp) <- names(ind1) <- colnames(x)
  est <- list(coefficients = coefTemp, a0 = a0, index = ind1, lambda0 = lambda0, 
              residuals = e1, sigma = sqrt(var(e1)), call = match.call(), options = list(post = post, 
                                                                                         intercept = intercept, control = control))
  class(est) <- c("rlassologit")
  return(est)
}


############################################################################################################################### 

################# Methods for logistic Lasso

#' Methods for S3 object \code{rlassologit}
#'
#' Objects of class \code{rlassologit} are constructed by \code{rlassologit} or \code{rlassologit.fit}.
#' \code{print.rlassologit} prints and displays some information about fitted \code{rlassologit} objects.
#' \code{summary.rlassologit} summarizes information of a fitted \code{rlassologit} object.
#' \code{predict.rlassologit} predicts values based on a \code{rlassologit} object.
#' \code{model.matrix.rlassologit} constructs the model matrix of a lasso object.
#' @param object an object of class \code{rlassologit}
#' @param x an object of class \code{rlassologit}
#' @param all logical, indicates if coefficients of all variables (TRUE) should be displayed or only the non-zero ones (FALSE)
#' @param digits significant digits in printout
#' @param type type of prediction required. The default ('response) is on the scale of the response variable; the alternative 'link' is on the scale of the linear predictors.
#' @param newdata new data set for prediction
#' @param ... arguments passed to the print function and other methods
#' @keywords methods rlassologit
#' @rdname methods.rlassologit
#' @aliases methods.rlassologit print.rlassologit summary.rlassologit predict.rlassologit model.matrix.rlassologit
#' @export

predict.rlassologit <- function(object, newdata = NULL, type = "response", 
                                ...) {
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
  } else {
    varcoef <- names(object$coefficients)
    if (all(is.element(varcoef, colnames(newdata)))) {
      X <- as.matrix(newdata[, varcoef])
    } else {
      X <- as.matrix(newdata)
      stopifnot(ncol(X) == length(object$coefficients))
    }
  }
  n <- dim(X)[1]  #length(object$residuals)
  beta <- object$coefficients
  if (object$options[["intercept"]]) {
    yp <- object$a0 + X %*% as.vector(beta)
    if (dim(X)[2] == 0) 
      yp <- rep(object$a0, n)
    if (type == "response") 
      yhat <- 1/(1 + exp(-yp))
    if (type == "link") 
      yhat <- yp
  }
  if (!object$options[["intercept"]]) {
    yp <- X %*% as.vector(beta)
    if (dim(X)[2] == 0) 
      yp <- rep(0, n)
    yhat <- 1/(1 + exp(-yp))
    if (type == "response") 
      yhat <- 1/(1 + exp(-yp))
    if (type == "link") 
      yhat <- yp
  }
  return(yhat)
}

#' @rdname methods.rlassologit
#' @export

model.matrix.rlassologit <- function(object, ...) {
  
  # falls formula
  if (is.call(object$call[[2]])) {
    # falls kein Datensatz uebergeben
    if (is.null(object$call$data)) {
      X <- model.frame(object$call[[2]])
      mt <- attr(X, "terms")
      attr(mt, "intercept") <- 0
      mm <- model.matrix(mt, X)
      # falls Datensatz uebergeben
    } else {
      dataev <- eval(object$call$data)
      mm <- as.matrix(dataev[, names(object$coefficients)])
    }
  } else {
    # falls default
    mm <- eval(object$call[[2]])
  }
  return(mm)
}


#' @rdname methods.rlassologit
#' @export

print.rlassologit <- function(x, all = TRUE, digits = max(3L, getOption("digits") - 
                                                            3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    coeffs <- coef(x)
    if (x$options$intercept) {
      coeffs <- c(x$a0, coeffs)
      names(coeffs)[1] <- "(Intercept)"
      index <- cbind(1, x$index)
    }
    if (all) {
      cat("Coefficients:\n")
      print.default(format(coeffs, digits = digits), print.gap = 2L, 
                    quote = FALSE)
    } else {
      print.default(format(coeffs[index], digits = digits), print.gap = 2L, 
                    quote = FALSE)
    }
  } else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

#' @rdname methods.rlassologit
#' @export

summary.rlassologit <- function(object, all = TRUE, digits = max(3L, getOption("digits") - 
                                                                   3L), ...) {
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  cat("\nPost-Lasso Estimation: ", paste(deparse(object$options$post), 
                                         sep = "\n", collapse = "\n"), "\n", sep = " ")
  coefs <- object$coefficients
  p <- length(coefs)
  num.selected <- sum(abs(object$coefficients) > 0)
  cat("\nTotal number of variables:", p)
  cat("\nNumber of selected variables:", num.selected, "\n", sep = " ")
  cat("\n")
  if (all) {
    coefm <- matrix(NA, length(coefs), 1)
    coefm[, 1] <- coefs
    colnames(coefm) <- "Estimate"
    rownames(coefm) <- names(coefs)
    printCoefmat(coefm, digits = digits, na.print = "NA")
  } else {
    coefs <- coefs[abs(coefs) > 0]
    coefm <- matrix(NA, length(coefs), 1)
    coefm[, 1] <- coefs
    colnames(coefm) <- "Estimate"
    rownames(coefm) <- names(coefs)
    printCoefmat(coefm, digits = digits, na.print = "NA")  #,  P.values=TRUE, has.Pvalue=TRUE)
  }
  cat("\n")
  invisible(object)
}
