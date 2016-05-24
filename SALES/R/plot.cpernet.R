###############################################################
## This function is adapted/modified based on the plot
#    function from the glmnet package:
## Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via
#    Coordinate Descent.
## Journal of Statistical Software, 33(1), 1-22.
## URL http://www.jstatsoft.org/v33/i01/.
###############################################################

plot.cpernet <- function(x, xvar = c("norm", "lambda"), 
  color = FALSE, label = FALSE, ...) {
  xvar <- match.arg(xvar)
  lambda <- x$lambda
  par(mfrow = c(2, 1))
  ##beta should be in 'dgCMatrix' format
  beta <- x$beta
  df <- x$df.beta
  if (all(df == 0)) stop("All mean coefficients are zero")
  which <- nonzero(beta)
  beta <- as.matrix(beta[which, , drop = FALSE])
  switch(xvar, norm = {
    index <- apply(abs(beta), 2, sum)
    iname <- "L1 Norm"
    }, lambda = {
         index <- log(lambda)
         iname <- "Log Lambda"})
  xlab <- iname
  ylab <- "Coefficients"
  dotlist <- list(...)
  type <- dotlist$type
  if (is.null(type)) {
    if (color == FALSE) 
      matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, 
        type = "l", pch = 500, col = gray.colors(12, start = 0.05, 
        end = 0.7, gamma = 2.2), ...) else matplot(index, t(beta),
        lty = 1, xlab = xlab, ylab = ylab, type = "l", pch = 500, ...)
  } else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, ...)
  atdf <- pretty(index)
  prettydf <- trunc(approx(x = index, y = df, xout = atdf, rule = 2)$y)
  axis(3, at = atdf, labels = prettydf, cex.axis = 1, tcl = NA)
  if (label) {
    nnz <- length(which)
    xpos <- max(index)
    pos <- 4
    if (xvar == "lambda") {
      xpos <- min(index)
      pos <- 2
    }
    xpos <- rep(xpos, nnz)
    ypos <- beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
  ##theta should be in 'dgCMatrix' format
  theta <- x$theta
  df <- x$df.theta
  if (all(df == 0)) stop("All scale coefficients are zero")
  which <- nonzero(theta)
  theta <- as.matrix(theta[which, , drop = FALSE])
  switch(xvar, norm = {
    index <- apply(abs(theta), 2, sum)
    iname <- "L1 Norm"
    }, lambda = {
         index <- log(lambda)
         iname <- "Log Lambda"})
  xlab <- iname
  ylab <- "Coefficients"
  dotlist <- list(...)
  type <- dotlist$type
  if (is.null(type)) {
    if (color == FALSE) 
      matplot(index, t(theta), lty = 1, xlab = xlab, ylab = ylab, 
        type = "l", pch = 500, col = gray.colors(12, start = 0.05, 
        end = 0.7, gamma = 2.2), ...) else matplot(index, t(theta),
        lty = 1, xlab = xlab, ylab = ylab, type = "l", pch = 500, ...)
  } else matplot(index, t(theta), lty = 1, xlab = xlab, ylab = ylab, ...)
  atdf <- pretty(index)
  prettydf <- trunc(approx(x = index, y = df, xout = atdf, rule = 2)$y)
  axis(3, at = atdf, labels = prettydf, cex.axis = 1, tcl = NA)
  if (label) {
    nnz <- length(which)
    xpos <- max(index)
    pos <- 4
    if (xvar == "lambda") {
      xpos <- min(index)
      pos <- 2
    }
    xpos <- rep(xpos, nnz)
    ypos <- beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
} 
