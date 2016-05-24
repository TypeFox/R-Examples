#' rigorous Lasso for Linear Models: Inference
#'
#' Estimation and inference of (low-dimensional) target coefficients in a high-dimensional linear model.
#'
#' The functions estimates (low-dimensional) target coefficients in a high-dimensional linear model.
#' An application is e.g. estimation of a treatment effect \eqn{\alpha_0} in a
#' setting of high-dimensional controls. The user can choose between the so-called post-double-selection method and partialling-out.
#' The idea of the double selection method is to select variables by Lasso regression of
#' the outcome variable on the control variables and the treatment variable on
#' the control variables. The final estimation is done by a regression of the
#' outcome on the treatment effect and the union of the selected variables in
#' the first two steps. In partialling-out first the effect of the regressors on the outcome and the treatment variable is taken out by Lasso and then a regression of the residuals is conducted. The resulting estimator for \eqn{\alpha_0} is normal
#' distributed which allows inference on the treatment effect. It presents a wrap function for \code{rlassoEffect} 
#' which does inference for a single variable.
#'
#' @param x matrix of regressor variables serving as controls and potential
#' treatments. For \code{rlassoEffect} it contains only controls, for \code{rlassoEffects} both controls and potential treatments. For \code{rlassoEffects} it must have at least two columns.
#' @param y outcome variable (vector or matrix)
#' @param index vector of integers, logicals or variables names indicating the position (column) of
#' variables (integer case), logical vector of length of the variables (TRUE or FALSE) or the variable names of \code{x} which should be used for inference / as treatment variables.
#' @param method method for inference, either 'partialling out' (default) or 'double selection'. 
#' @param I3 For the 'double selection'-method the logical vector \code{I3} has same length as the number of variables in \code{x};
#' indicates if variables (TRUE) should be included in any case to the model and they are exempt from selection. These variables should not be included in the \code{index}; hence the intersection with \code{index} must be the empty set.
#' In the case of partialling out it is ignored.
#' @param post logical, if post Lasso is conducted with default \code{TRUE}.
#' @param \dots parameters passed to the \code{rlasso} function.
#' @return The function returns an object of class \code{rlassoEffects} with the following entries: \item{coefficients}{vector with estimated
#' values of the coefficients for each selected variable} \item{se}{standard error (vector)}
#' \item{t}{t-statistic} \item{pval}{p-value} \item{samplesize}{sample size of the data set} \item{index}{index of the variables for which inference is performed}
#' @references A. Belloni, V. Chernozhukov, C. Hansen (2014). Inference on
#' treatment effects after selection among high-dimensional controls. The
#' Review of Economic Studies 81(2), 608-650.
#' @keywords Estimation Inference Treatment effect High-dimensional controls
#' @export
#' @rdname rlassoEffects
#' @examples
#' library(hdm)
#' ## DGP
#' n <- 250
#' p <- 100
#' px <- 10
#' X <- matrix(rnorm(n*p), ncol=p)
#' beta <- c(rep(2,px), rep(0,p-px))
#' intercept <- 1
#' y <- intercept + X %*% beta + rnorm(n)
#' ## fit rlassoEffects object with inference on three variables
#' rlassoEffects.reg <- rlassoEffects(x=X, y=y, index=c(1,7,20))
#' ## methods
#' summary(rlassoEffects.reg)
#' confint(rlassoEffects.reg, level=0.9)
rlassoEffects <- function(x, y, index = c(1:ncol(x)), method = "partialling out", 
                          I3 = NULL, post = TRUE, ...) {
  
  checkmate::checkChoice(method, c("partialling out", "double selection"))
  
  if (is.logical(index)) {
    k <- p1 <- sum(index)
  } else {
    k <- p1 <- length(index)
  }
  n <- dim(x)[1]
  x <- as.matrix(x)
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
  if (method == "double selection") {
    # check validity of I3
    I3ind <- which(I3 == T)
    if (length(intersect(index, I3ind) != 0)) 
      stop("I3 and index must not overlap!")
  }
  
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
    lasso.regs[[i]] <- try(col <- rlassoEffect(Xt, y, d, method = method, 
                                               I3 = I3m, post = post, ...))
    if (class(lasso.regs[[i]]) == "try-error") {
      next
    } else {
      coefficients[i] <- col$alpha
      se[i] <- col$se
      t[i] <- col$t
      pval[i] <- col$pval
      reside[, i] <- col$residuals$epsilon
      residv[, i] <- col$residuals$v
    }
  }
  residuals <- list(e = reside, v = residv)
  res <- list(coefficients = coefficients, se = se, t = t, pval = pval, 
              lasso.regs = lasso.regs, index = index, call = match.call(), samplesize = n, 
              residuals = residuals)
  class(res) <- "rlassoEffects"
  return(res)
}

#' @rdname rlassoEffects
#' @param d variable for which inference is conducted (treatment variable)
#' @export
rlassoEffect <- function(x, y, d, method = "double selection", I3 = NULL, 
                         post = TRUE, ...) {
  d <- as.matrix(d, ncol = 1)
  y <- as.matrix(y, ncol = 1)
  kx <- dim(x)[2]
  n <- dim(x)[1]
  if (is.null(colnames(d))) 
    colnames(d) <- "d1"
  if (is.null(colnames(x)) & !is.null(x)) 
    colnames(x) <- paste("x", 1:kx, sep = "")
  if (method == "double selection") {
    I1 <- rlasso(d ~ x, post = post, ...)$index
    I2 <- rlasso(y ~ x, post = post, ...)$index
    
    
    if (is.logical(I3)) {
      I <- I1 + I2 + I3
      I <- as.logical(I)
    } else {
      I <- I1 + I2
      I <- as.logical(I)
    }
    if (sum(I) == 0) {
      I <- NULL
    }
    x <- cbind(d, x[, I, drop = FALSE])
    reg1 <- lm(y ~ x)
    alpha <- coef(reg1)[2]
    xi <- reg1$residuals * sqrt(n/(n - sum(I) - 1))
    if (is.null(I)) {
      reg2 <- lm(d ~ 1)
    }
    if (!is.null(I)) {
      reg2 <- lm(d ~ x[, -1, drop = FALSE])
    }
    v <- reg2$residuals
    var <- 1/n * 1/mean(v^2) * mean(v^2 * xi^2) * 1/mean(v^2)
    se <- sqrt(var)
    tval <- alpha/sqrt(var)
    pval <- 2 * pnorm(-abs(tval))
    if (is.null(I)) {
      no.selected <- 1
    } else {
      no.selected <- 0
    }
    res <- list(epsilon = xi, v = v)
    # results <- list(alpha=unname(alpha), se=drop(se), t=unname(tval),
    # pval=unname(pval), no.selected=no.selected,
    # coefficients=unname(alpha), coefficient=unname(alpha),
    # coefficients.reg=coef(reg1), residuals=res, call=match.call(),
    # samplesize=n)
    results <- list(alpha = alpha, se = drop(se), t = tval, pval = pval, 
                    no.selected = no.selected, coefficients = alpha, coefficient = alpha, 
                    coefficients.reg = coef(reg1), residuals = res, call = match.call(), 
                    samplesize = n)
  }
  
  if (method == "partialling out") {
    yr <- rlasso(y ~ x, post = post, ...)$residuals
    dr <- rlasso(d ~ x, post = post, ...)$residuals
    reg1 <- lm(yr ~ dr)
    alpha <- coef(reg1)[2]
    var <- vcov(reg1)[2, 2]
    se <- sqrt(var)
    tval <- alpha/sqrt(var)
    pval <- 2 * pnorm(-abs(tval))
    res <- list(epsilon = reg1$residuals, v = dr)
    results <- list(alpha = unname(alpha), se = drop(se), t = unname(tval), 
                    pval = unname(pval), coefficients = unname(alpha), coefficient = unname(alpha), 
                    coefficients.reg = coef(reg1), residuals = res, call = match.call(), 
                    samplesize = n)
  }
  class(results) <- "rlassoEffects"
  return(results)
}

################################################################################################################################### Methods for rlassoEffects

#' Methods for S3 object \code{rlassoEffects}
#'
#' Objects of class \code{rlassoEffects} are constructed by  \code{rlassoEffects}.
#' \code{print.rlassoEffects} prints and displays some information about fitted \code{rlassoEffect} objects.
#' summary.rlassoEffects summarizes information of a fitted \code{rlassoEffect} object and is described at \code{\link{summary.rlassoEffects}}.
#' \code{confint.rlassoEffects} extracts the confidence intervals.
#' \code{plot.rlassoEffects} plots the estimates with confidence intervals.
#'
#' @param x an object of class \code{rlassoEffects}
#' @param digits significant digits in printout
#' @param ... arguments passed to the print function and other methods.
#' @keywords methods rlassoEffects
#' @rdname methods.rlassoEffects
#' @aliases methods.rlassoEffects print.rlassoEffects confint.rlassoEffects plot.rlassoEffects
#' @export

print.rlassoEffects <- function(x, digits = max(3L, getOption("digits") - 
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

#' @rdname methods.rlassoEffects
#' @param object an object of class \code{rlassoEffects}
#' @param parm a specification of which parameters are to be given confidence intervals among the variables for which inference was done, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level required
#' @param joint logical, if \code{TRUE} joint confidence intervals are calculated.
#' @export

confint.rlassoEffects <- function(object, parm, level = 0.95, joint = FALSE, 
                                  ...) {
  B <- 500  # number of bootstrap repitions
  n <- object$samplesize
  k <- p1 <- length(object$coefficients)
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- pnames else if (is.numeric(parm)) 
      parm <- pnames[parm]
  if (!joint) {
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    # fac <- qt(a, n-k)
    fac <- qnorm(a)
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

#' @rdname methods.rlassoEffects
#' @export
#' @param main an overall title for the plot
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param xlim vector of length two giving lower and upper bound of x axis
plot.rlassoEffects <- function(x, main = "", xlab = "coef", ylab = "", 
                               xlim = NULL, ...) {
  
  # generate ordered KI-matrix
  coefmatrix <- cbind(summary(x)$coef, confint(x))[, c(1, 5, 6)]
  if (is.null(dim(coefmatrix))) {
    vec <- coefmatrix
    coefmatrix <- matrix(vec, ncol = 3)
    colnames(coefmatrix) <- names(vec)
  }
  
  rownames(coefmatrix) <- names(x$coefficients)
  coefmatrix <- as.data.frame(coefmatrix)
  coefmatrix <- cbind(rownames(coefmatrix), coefmatrix)
  colnames(coefmatrix) <- c("names", "coef", "lower", "upper")
  coefmatrix <- coefmatrix[order(abs(coefmatrix[, 2])), ]
  
  col <- "#000099"
  # scale
  if (missing(xlim)) {
    low <- min(coefmatrix[, -1])
    up <- max(coefmatrix[, -1])
  } else {
    low <- xlim[1]
    up <- xlim[2]
  }
  # generate points
  plotobject <- ggplot2::ggplot(coefmatrix, ggplot2::aes(y = coef, x = factor(names, 
                                                                              levels = names))) + ggplot2::geom_point(colour = col, size = 1.75) + 
    ggplot2::geom_hline(h = 0, colour = col, width = 0.1)
  
  # generate errorbars (KIs)
  plotobject <- plotobject + ggplot2::geom_errorbar(ymin = coefmatrix$lower, 
                                                    ymax = coefmatrix$upper, colour = col, width = 0.4, size = 0.2)
  
  # further graphic parameter
  plotobject <- plotobject + ggplot2::ggtitle(main) + ggplot2::ylim(low, 
                                                                    up) + ggplot2::xlab(ylab) + ggplot2::ylab(xlab)
  
  
  ## invert x and y axis
  #plotobject <- plotobject + ggplot2::coord_flip()
  
  # layout
  plotobject <- plotobject + ggplot2::theme_bw() + ggplot2::geom_blank() + 
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank())
  # plot
  plotobject
}


################ Methods: summary

#' Summarizing rlassoEffects fits
#' 
#' Summary method for class \code{rlassoEffects}
#' 
#' Summary of objects of class \code{rlassoEffects}
#' 
#' @param object an object of class \code{rlassoEffects}, usually a result of a call to \code{rlassoEffects}
#' @param ... further arguments passed to or from other methods.
#' @rdname summary.rlassoEffects
#' @export 
summary.rlassoEffects <- function(object, ...) {
  ans <- NULL
  k <- length(object$coefficients)
  table <- matrix(NA, ncol = 4, nrow = k)
  rownames(table) <- names(object$coefficient)
  colnames(table) <- c("Estimate.", "Std. Error", "t value", "Pr(>|t|)")
  table[, 1] <- object$coefficients
  table[, 2] <- object$se
  table[, 3] <- object$t
  table[, 4] <- object$pval
  ans$coefficients <- table
  ans$object <- object
  class(ans) <- "summary.rlassoEffects"
  return(ans)
}


#' @param x an object of class \code{summary.rlassoEffects}, usually a result of a call or \code{summary.rlassoEffects}
#' @param digits the number of significant digits to use when printing.
#' @method print summary.rlassoEffects
#' @rdname summary.rlassoEffects
#' @export
print.summary.rlassoEffects <- function(x, digits = max(3L, getOption("digits") - 
                                                          3L), ...) {
  if (length(coef(x$object))) {
    k <- dim(x$coefficients)[1]
    table <- x$coefficients
    print("Estimates and significance testing of the effect of target variables")
    printCoefmat(table, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}