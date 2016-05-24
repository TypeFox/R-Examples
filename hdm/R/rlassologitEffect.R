#' rigorous Lasso for Logistic Models: Inference
#'
#' The function estimates (low-dimensional) target coefficients in a high-dimensional logistic model.
#'
#' The functions estimates (low-dimensional) target coefficients in a high-dimensional logistic model.
#' An application is e.g. estimation of a treatment effect \eqn{\alpha_0} in a
#' setting of high-dimensional controls. The function is a wrap function for \code{rlassologitEffect} which does inference for only one variable (d).
#'
#' @param x matrix of regressor variables serving as controls and potential
#' treatments.  For \code{rlassologitEffect} it contains only controls, for \code{rlassologitEffects} both controls and potential treatments. For  \code{rlassologitEffects} it must have at least two columns.
#' @param y outcome variable
#' @param index vector of integers, logical or names indicating the position (column) or name of
#' variables of x which should be used as treatment variables.
#' @param I3 logical vector with same length as the number of controls;
#' indicates if variables (TRUE) should be included in any case.
#' @param \dots additional parameters
#' @return The function returns an object of class \code{rlassologitEffects} with the following entries: \item{coefficients}{estimated
#' value of the coefficients} \item{se}{standard errors}
#' \item{t}{t-statistics} \item{pval}{p-values} \item{samplesize}{sample size of the data set} \item{I}{index of variables of the union of the lasso regressions}
#' @references A. Belloni, V. Chernozhukov, Y. Wei (2013). Honest confidence regions for a regression parameter in logistic regression with a loarge number of controls.
#' cemmap working paper CWP67/13.
#' @keywords Estimation Inference Logistic Lasso
#' @export
#' @rdname rlassologitEffects
#' @export
rlassologitEffects <- function(x, y, index = c(1:ncol(x)), I3 = NULL, ...) {
  if (is.logical(index)) {
    k <- p1 <- sum(index)
  } else {
    k <- p1 <- length(index)
  }
  n <- dim(x)[1]
  # preprocessing index numerischer Vektor
  if (is.numeric(index)) {
    index <- as.integer(index)
    stopifnot(all(index <= ncol(x)) && length(index) <= ncol(x))
  } else {
    # logical Vektor
    if (is.logical(index)) {
      stopifnot(length(index) == ncol(x))
      index <- which(index == T)
    } else {
      # character Vektor
      if (is.character(index)) {
        stopifnot(all(is.element(index, colnames(x))))
        index <- which(is.element(colnames(x), index))
      } else {
        stop("argument index has an invalid type")
      }
    }
  }
  
  # check validity of I3
  I3ind <- which(I3 == T)
  if (length(intersect(index, I3ind) != 0)) 
    stop("I3 and index must not overlap!")
  if (is.null(colnames(x))) 
    colnames(x) <- paste("V", 1:dim(x)[2], sep = "")
  coefficients <- as.vector(rep(NA_real_, k))
  se <- rep(NA_real_, k)
  t <- rep(NA_real_, k)
  pval <- rep(NA_real_, k)
  lasso.regs <- vector("list", k)
  reside <- matrix(NA, nrow = n, ncol = p1)
  residv <- matrix(NA, nrow = n, ncol = p1)
  names(coefficients) <- names(se) <- names(t) <- names(pval) <- names(lasso.regs) <- colnames(reside) <- colnames(residv) <- colnames(x)[index]
  
  for (i in 1:k) {
    d <- x[, index[i], drop = FALSE]
    Xt <- x[, -index[i], drop = FALSE]
    I3m <- I3[-index[i]]
    lasso.regs[[i]] <- try(col <- rlassologitEffect(Xt, y, d, I3 = I3m))
    if (class(lasso.regs[[i]]) == "try-error") {
      next
    } else {
      coefficients[i] <- col$alpha
      se[i] <- col$se
      t[i] <- col$t
      pval[i] <- col$pval
      # reside[,i] <- col$residuals$epsilon residv[,i] <- col$residuals$v
    }
  }
  residuals <- list(e = reside, v = residv)
  res <- list(coefficients = coefficients, se = se, t = t, pval = pval, 
              lasso.regs = lasso.regs, index = I, call = match.call(), samplesize = n, 
              residuals = residuals)
  class(res) <- "rlassologitEffects"
  return(res)
}


#' @rdname rlassologitEffects
#' @param d variable for which inference is conducted (treatment variable)
#' @export
#'
rlassologitEffect <- function(x, y, d, I3 = NULL) {
  d <- as.matrix(d, ncol = 1)
  y <- as.matrix(y, ncol = 1)
  kx <- p <- dim(x)[2]
  n <- dim(x)[1]
  if (is.null(colnames(d))) 
    colnames(d) <- "d1"
  if (is.null(colnames(x)) & !is.null(x)) 
    colnames(x) <- paste("x", 1:kx, sep = "")
  # Step 1
  la1 <- 1.1/2 * sqrt(n) * qnorm(1 - 0.05/(max(n, (p + 1) * log(n))))
  dx <- cbind(d, x)
  l1 <- rlassologit(y ~ dx, post = TRUE, intercept = TRUE, penalty = list(lambda.start = la1))
  t <- predict(l1, type = "link", newdata = dx)
  sigma2 <- exp(t)/(1 + exp(t))^2
  w <- sigma2  #exp(t)/(1+exp(t))^2
  f <- w/sigma2
  I1 <- l1$index[-1]
  # Step 2
  la2 <- rep(2.2 * sqrt(n) * qnorm(1 - 0.05/(max(n, p * log(n)))), p)
  xf <- x * as.vector(f)
  df <- d * f
  l2 <- rlasso(xf, df, post = TRUE, intercept = TRUE, penalty = list(homoscedastic = "none", 
                                                                     lambda.start = la2, c = 1.1, gamma = 0.1))
  I2 <- l2$index
  z <- l2$residual/sqrt(sigma2)
  # Step 3
  if (is.logical(I3)) {
    I <- I1 + I2 + I3
    I <- as.logical(I)
  } else {
    I <- I1 + I2
    I <- as.logical(I)
  }
  xselect <- x[, I]
  p3 <- dim(xselect)[2]
  la3 <- 1.1/2 * sqrt(n) * qnorm(1 - 0.05/(max(n, (p3 + 1) * log(n))))
  l3 <- rlassologit(cbind(d, xselect), y, post = TRUE, normalize = TRUE, 
                    intercept = TRUE, penalty = list(lambda.start = la3))
  alpha <- l3$coefficients[1]
  t3 <- predict(l3, type = "link", newdata = cbind(d, xselect))
  G3 <- exp(t3)/(1 + exp(t3))
  w3 <- G3 * (1 - G3)
  S21 <- 1/mean(w3 * d * z)^2 * mean((y - G3)^2 * z^2)
  xtilde <- x[, l3$index[-1]]
  p2 <- sum(l3$index) + 1
  b <- cbind(d, xtilde)
  p4 <- dim(b)[2]
  A <- matrix(0, ncol = p4, nrow = p4)
  for (i in 1:n) {
    A <- A + w3[i] * outer(b[i, ], b[i, ])
  }
  S22 <- solve(1/n * A)[1, 1]
  S2 <- max(S21, S22)
  se <- sqrt(S2/n)
  tval <- alpha/se
  pval <- 2 * pnorm(-abs(tval))
  if (is.null(I)) {
    no.selected <- 1
  } else {
    no.selected <- 0
  }
  # return(list(alpha=unname(alpha), se=drop(se), t=unname(tval),
  # pval=unname(pval), coefficients=coef(l3), residuals=l3$residuals))
  results <- list(alpha = alpha, se = drop(se), t = tval, pval = pval, 
                  no.selected = no.selected, coefficients = alpha, coefficient = alpha, 
                  residuals = l3$residuals, call = match.call(), samplesize = n)
  return(results)
}



################# Methods for rlassologitEffects

#' Methods for S3 object \code{rlassologitEffects}
#'
#' Objects of class \code{rlassologitEffects} are construced by \code{rlassologitEffects} or \code{rlassologitEffect}. 
#' \code{print.rlassologitEffects} prints and displays some information about fitted \code{rlassologitEffect} objects.
#' \code{summary.rlassologitEffects} summarizes information of a fitted \code{rlassologitEffects} object.
#' \code{confint.rlassologitEffects} extracts the confidence intervals.
#' @param object an object of class \code{rlassologitEffects}
#' @param x an object of class \code{rlassologitEffects}
#' @param digits number of significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @keywords methods rlassologitEffects
#' @rdname methods.rlassologitEffects
#' @aliases methods.rlassologitEffects print.rlassologitEffects summary.rlassologitEffects confint.rlassologitEffects
#' @export

print.rlassologitEffects <- function(x, digits = max(3L, getOption("digits") - 
                                                       3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L, 
                  quote = FALSE)
  } else cat("No coefficients\n")
  cat("\n")
  invisible(coef(x))
}


#' @rdname methods.rlassologitEffects
#' @export

summary.rlassologitEffects <- function(object, digits = max(3L, getOption("digits") - 
                                                              3L), ...) {
  if (length(coef(object))) {
    k <- length(object$coefficients)
    table <- matrix(NA, ncol = 4, nrow = k)
    rownames(table) <- names(object$coefficient)
    colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
    table[, 1] <- object$coefficient
    table[, 2] <- object$se
    table[, 3] <- object$t
    table[, 4] <- object$pval
    print("Estimates and significance testing of the effect of target variables")
    printCoefmat(table, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}

#' @rdname methods.rlassologitEffects
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level required.
#' @param joint logical, if joint confidence intervals should be clalculated
#' @export

confint.rlassologitEffects <- function(object, parm, level = 0.95, joint = FALSE, 
                                       ...) {
  B <- 500  # number of bootstrap repitions
  n <- object$samplesize
  k <- p1 <- length(object$coefficient)
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- pnames else if (is.numeric(parm)) 
      parm <- pnames[parm]
  if (!joint) {
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    fac <- qt(a, n - k)
    pct <- format.perc(a, 3)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                               pct))
    ses <- object$se[parm]
    ci[] <- cf[parm] + ses %o% fac
  }
  
  if (joint) {
    phi <- object$residuals$e * object$residuals$v
    m <- 1/sqrt(colMeans(phi^2))
    phi <- t(t(phi)/m)
    sigma <- sqrt(colMeans(phi^2))
    sim <- vector("numeric", length = B)
    for (i in 1:B) {
      xi <- rnorm(n)
      phi_temp <- phi * xi
      Nstar <- 1/sqrt(n) * colSums(phi_temp)
      sim[i] <- max(abs(Nstar))
    }
    a <- (1 - level)/2
    ab <- c(a, 1 - a)
    pct <- format.perc(ab, 3)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                                               pct))
    hatc <- quantile(sim, probs = 1 - a)
    ci[, 1] <- cf[parm] - hatc * 1/sqrt(n) * sigma
    ci[, 2] <- cf[parm] + hatc * 1/sqrt(n) * sigma
  }
  return(ci)
}





