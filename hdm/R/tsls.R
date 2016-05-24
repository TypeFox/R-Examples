#' Two-Stage Least Squares Estimation (TSLS)
#'
#' The function does Two-Stage Least Squares Estimation (TSLS).
#'
#' The function computes tsls estimate (coefficients) and variance-covariance-matrix assuming homoskedasticity
#' for outcome variable \code{y} where \code{d} are endogenous variables in structural equation, \code{x} are exogensous variables in
#' structural equation and z are instruments. It returns an object of class \code{tsls} for which the methods \code{print} and \code{summary} 
#' are provided.
#'
#' @param y outcome variable
#' @param d endogenous variables
#' @param x exogenous variables
#' @param z instruments
#' @param intercept logical, if intercept should be included
#' @param homoscedastic logical, if homoscedastic (\code{TRUE}, default) or heteroscedastic erros (\code{FALSE}) should be calculated.
#' @return The function returns a list with the following elements \item{coefficients}{coefficients}
#' \item{vcov}{variance-covariance matrix} \item{residuals}{outcome minus predicted values} \item{call}{function call} \item{samplesize}{sample size}
#' \item{se}{standard error}
#' @keywords Instrumental Variables
#' @keywords Endogeneity
#' @keywords 2SLS
#' @keywords TSLS
#' @export


tsls <- function(y, d, x, z, intercept=TRUE, homoscedastic=TRUE) {
  n <- length(y)
  
  d <- as.matrix(d)
  if (!is.null(x))  x <- as.matrix(x)
  z <- as.matrix(z)
  if (is.null(colnames(d))  & is.matrix(d)) colnames(d) <- paste("d", 1:ncol(d), sep="")
  if (is.null(colnames(x)) & !is.null(x) & is.matrix(x)) colnames(x) <- paste("x", 1:ncol(x), sep="")
  if (is.null(colnames(z)) & !is.null(z) & is.matrix(z)) colnames(z) <- paste("z", 1:ncol(z), sep="")
  
  if (intercept==TRUE && is.null(x)) {
    x <- as.matrix(rep(1,n))
    colnames(x) <- "Intercept"
  } else {
    if (intercept==TRUE) {
      x <- as.matrix(cbind(1,x))
      colnames(x)[1] <- "Intercept"
    }
  }
  
  a1 <- dim(d)[2]
  a2 <- dim(x)[2]
  if (is.null(x)) {
    a2 <- 0
  }
  if (is.vector(x)) {
    a2 <- 1
  }
  if (is.vector(d)) {
    a1 <- 1
  }
  k <- a1 + a2
  X <- cbind(d, x)
  Z <- cbind(z, x)

  Mxz <- t(X) %*% Z
  Mzz <- solve(t(Z) %*% Z)
  M <- solve(Mxz %*% Mzz %*% t(Mxz))

  b <- M %*% Mxz %*% Mzz %*% (t(Z) %*% y)
  #Dhat <- Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% X
  #b2 <- MASS::ginv(t(Dhat) %*% X) %*% (t(Dhat) %*% y)
  if (homoscedastic==TRUE) {
  e <- y - X %*% b
  VC1 <- as.numeric((t(e) %*% e/(n - k))) * M
  }
  if (homoscedastic==FALSE) {
    e <- y - X %*% b
    S <- 0
    for (i in 1:n) {
      S <- S + e[i]^2*(Z[i,]%*%t(Z[i,]))
    }
    S <- S*1/n
    VC1 <- n*M%*%(Mxz%*%Mzz%*%S%*%Mzz%*%t(Mxz))%*%M
  }
  res <- list(coefficients = b, vcov = VC1, se=sqrt(diag(VC1)), residuals = e, call=match.call(), samplesize=n)
  class(res) <- "tsls"
  return(res)
}

################# Methods for tsls

#' Methods for S3 object \code{tsls}
#'
#' Objects of class \code{tsls} are constructed by \code{tsls}. 
#' \code{print.tsls} prints and displays some information about fitted \code{tsls} objects.
#' \code{summary.tsls} summarizes information of a fitted \code{tsls} object.
#' @param x an object of class \code{tsls}
#' @param digits significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @keywords methods tsls
#' @rdname methods.tsls
#' @aliases methods.tsls print.tsls summary.tsls
#' @export

print.tsls <- function(x, digits = max(3L, getOption("digits") - 
                                                    3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    coeffs <- as.matrix(coef(x))
    colnames(coeffs) <- "Estimate"
    cat("Coefficients:\n")
    print.default(format(coeffs, digits = digits), print.gap = 2L, 
                  quote = FALSE)
  } else cat("No coefficients\n")
  cat("\n")
  invisible(coef(x))
}

#' @param object an object of class \code{tsls}
#' @rdname methods.tsls
#' @export

summary.tsls <- function(object, digits = max(3L, getOption("digits") - 
                                                           3L), ...) {
  if (length(coef(object))) {
    k <- length(object$coefficient)
    table <- matrix(NA, ncol = 4, nrow = k)
    rownames(table) <- dimnames(object$coefficients)[[1]] #names(object$coefficient)
    colnames(table) <- c("Estimate", "Std. Error", "t value", "p value")
    table[, 1] <- object$coefficient
    table[, 2] <- sqrt(diag(as.matrix(object$vcov)))
    table[, 3] <- table[, 1]/table[, 2]
    table[, 4] <- 2 * pnorm(-abs(table[, 3]))
    print("Estimates and Significance Testing from from tsls")
    printCoefmat(table, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}