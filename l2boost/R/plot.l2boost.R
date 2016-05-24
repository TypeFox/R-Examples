#---------------------------------------------------------------------
# l2boost plots
# ... currently allows for rho/coef paths and will also plot cv objects
# ... standarized rho paths are gradient-correlation paths
# ... standarized coefficient paths are coefficient paths for standardized x
# SHOULD WORK FOR ALL l2boost variants
#---------------------------------------------------------------------

#' @title Plotting for \code{\link{l2boost}} objects.
#'
#' @description 
#' plotting methods for \code{\link{l2boost}} objects (\code{\link{l2boost}} and \code{\link{cv.l2boost}}). 
#' 
#' By default, plotting an \code{\link{l2boost}} object produces a gradient-correlation vs iteration steps (m) plot.
#' Plotting a \code{\link{cv.l2boost}} object produces a cross-validation error plot, and prints the minimal CV MSE value
#' and optimal step opt.step to the R console.
#' 
#' Many generic arguments to \code{\link{plot}} are passed through the \code{\link{plot.l2boost}} function.
#' 
#' @details
#' Gradient-correlation plots are created by tracing out the boosting coefficient (rho) for each candidate
#' direction. The coefficient and gradient-correlation are equivalent under standard scaling (zero intercept with 
#' design matrix columns scaled to have mean=0 and variance=1).  
#' 
#' Unless explicitly set using \emph{col} argument, the plot function colors the gradient-correlation paths along each
#' direction by the following criteria: 
#' \itemize{
#' \item Red: indicates the coordinate direction has been selected in the boosting path at some step <= m. 
#' \item Blue: indicates the coordinate will be selected within the specified number of steps M (and switch to 
#' red upon selection).
#' \item Grey: indicates coordinates have not and will not be selected by the algorithm over all iterations. 
#' }
#' The colors are set using the \emph{l.crit} return value from the \code{\link{l2boost}} object.
#' 
#' @param x l2boost or cv.l2boost object
#' @param type which type of plot. \emph{rho} plots gradient-correlation, \emph{coef} regression (beta) 
#' coefficients vs the step number m along the x-axis
#' @param standardize Should we plot standardized gradient-correlation (default: TRUE)
#' @param active.set Vector of indices of the coordinates for highlighting with 
#' color=col (default: NULL shows all active coordinates)
#' @param xvar what measure do we plot on the x-axis? \emph{step} plots the step m, \emph{norm} plots the 
#' normalized distance (1-nu)^(m-1)
#' @param xlab specific x-axis label (NULL results in default value depending on xvar)
#' @param ylab specific y-axis label (NULL results in default value depending on type)
#' @param trim (default: TRUE)
#' @param clip Do we want to c
#' @param col Color to highlight active.set coordinates (NULL indicates default all active set at 
#' step M in blue, changes to red after selection
#' @param ylim Control plotted y-values (default: NULL for auto range)
#' @param xlim Control plotted x-values (default: NULL for auto domain )
#' @param ... other arguments passed to plot functions
#'
#' @return \code{NULL}
#'
#' @seealso \code{\link{l2boost}}, \code{\link{print.l2boost}}, \code{\link{predict.l2boost}} methods of l2boost 
#' and \code{\link{cv.l2boost}}
#'
#' @examples
#' #--------------------------------------------------------------------------
#' # Example: Diabetes 
#' #  
#' # See Efron B., Hastie T., Johnstone I., and Tibshirani R. 
#' # Least angle regression. Ann. Statist., 32:407-499, 2004.
#' data(diabetes, package = "l2boost")
#' 
#' l2.object <- l2boost(diabetes$x,diabetes$y, M=1000, nu=.01)
#'
#' # Plot the gradient-correlation, and regression beta coefficients as a function of
#' # boosting steps m
#' par(mfrow=c(2,2))
#' plot(l2.object)
#' abline(v=500, lty=2, col="grey")
#' plot(l2.object, type="coef")
#' abline(v=500, lty=2, col="grey")
#' 
#' # limit the plot to only the first 500 steps of the algorithm 
#' # (grey vertical line in previous plots).
#' plot(l2.object, xlim=c(0,500))
#' plot(l2.object, type="coef", xlim=c(0,500))
#' 
#' \dontrun{
#' #--------------------------------------------------------------------------
#' # Example: Plotting cross-validation objects
#' dta <- elasticNetSim(n=100)
#' # Set the boosting parameters
#' Mtarget = 1000
#' nuTarget = 1.e-2
#' 
#' cv.l2 <- cv.l2boost(dta$x,dta$y,M=Mtarget, nu=nuTarget, lambda=NULL)
#' 
#' # Show the CV MSE plot, with a marker at the "optimal iteration"
#' plot(cv.l2)
#' abline(v=cv.l2$opt.step, lty=2, col="grey")
#' 
#' # Show the l2boost object plots.
#' plot(cv.l2$fit)
#' abline(v=cv.l2$opt.step, lty=2, col="grey")
#'  
#' plot(cv.l2$fit, type="coef")
#' abline(v=cv.l2$opt.step, lty=2, col="grey")
#' 
#' # Create a color vector of length p=40 (from elasticNetSim defaults)
#' clr <- rep("black", 40)
#' # Set coordinates in the boosting path to color red.
#' clr[unique(cv.l2$fit$l.crit)] = "red"
#' 
#' # Show the "optimal" coefficient values, 
#' # red points are selected in boosting algorithm.
#' plot(coef(cv.l2$fit, m=cv.l2$opt.step), col=clr, ylab=expression(beta))
#' }
#' @method plot l2boost
#' @S3method plot l2boost
#' 
plot.l2boost <- function(x, 
                         type = c("rho", "coef"),
                         standardize = TRUE, active.set=NULL,
                         xvar = c("step", "norm"),
                         xlab = NULL, ylab = NULL,
                         trim = TRUE, clip=NULL, col=NULL,ylim=NULL, xlim=NULL,...) {
  # preliminary checks
  if (class(x)[1] != "l2boost") stop("This function only works for xs of class `l2boost'")
  type <- match.arg(type)
  xvar <- match.arg(xvar)
  # ----------------------------------------------------------------------------
  # cv plots
  # ----------------------------------------------------------------------------
  if (class(x)[2] == "cv") {
    # convert mse into matrix format more conducive for plotting/printing
    mse <- x$mse.list
    K <- length(mse)
    M <- max(sapply(1:K, function(k){length(mse[[k]])}), na.rm = TRUE)
    cv.all <-  matrix(NA, nrow = M, ncol = K)
    x.range <- rep(NA, K)
    for (k in 1:K) {
      cv.all[1:length(mse[[k]]), k] <- mse[[k]]
      x.range[k] <- length(mse[[k]])
    }
    step.size <- c(1, max(x.range, na.rm = TRUE))
    cv <- apply(cv.all, 1, mean, na.rm = TRUE)
    cv.error <- sqrt(apply(cv.all, 1, VAR)/K)
    #trim the y-axis
    if (trim) {
      if (K <= 5) {
        y.range <- quantile(c(cv.all, cv, cv + cv.error, cv - cv.error), c(0, .95), na.rm = TRUE)
      }
      else {
        y.range <- quantile(c(cv, cv + cv.error, cv - cv.error), c(0, .95), na.rm = TRUE)
      }
    }
    else {
      if (K <= 5) {
        y.range <- range(c(cv.all, cv, cv + cv.error, cv - cv.error), na.rm = TRUE)
      }
      else {
        y.range <- range(c(cv, cv + cv.error, cv - cv.error), na.rm = TRUE)
      }
    }
    if (is.null(xlab)) xlab <- "steps"
    if (is.null(ylab)) ylab <- "cross-validated MSE"
    matplot(1:M, cv.all, type = c("l", "n")[1 + 1 * (K > 5)], lty = 3, col = 2, lwd = 0.05,
            xlim = range(step.size),
            ylim = y.range,
            xlab = xlab, ylab = ylab, ...=...)
    lines(1:M, cv, lty = 1, lwd = 5, col = 2)
    error.bars(1:M, cv + cv.error, cv - cv.error, width = 0.0025, col = "gray")
    cat("minimum cross-validated MSE equals", round(x$mse, 4), "for step size", 
        x$opt.step -1, "\n")
  }
  else {
    y <- x$y
    Fm.path <- x$Fm.path
    rhom.path <- x$rhom.path
    M <- length(rhom.path)
    p <- length(x$betam)
    l.crit <- x$l.crit
    
    # determine what goes on the x-axis: (i) step (ii) norm
    if (xvar == "step") {
      xval <- 1:M
      if (is.null(xlab)) xlab <- "step"
    }
    else {
      b.m.path <- predict.l2boost(x, type = "coef")$coef.path
      xval <- sapply(1:M, function(m) {sum(abs(b.m.path[[m]]), na.rm = TRUE)})
      if (is.null(xlab)) xlab <- "l1-norm"
    }
    
    # ----------------------------------------------------------------------------
    # rho path plots
    # ----------------------------------------------------------------------------
    if (type == "rho") {
      if (is.null(ylab)) ylab <- "gradient"
      if (standardize) {
        rhom.path <- lapply(1:M, function(m){rhom.path[[m]]/sqrt(sum((y - Fm.path[[m]])^2, na.rm = TRUE))})
        if (is.null(ylab)) ylab <- "gradient-correlation"
      }
      path <- rhom.path
    }
    # ----------------------------------------------------------------------------
    # coef path plots
    # ----------------------------------------------------------------------------
    else if (type == "coef") {
      if (standardize) {
        if (is.null(ylab)) ylab <- "standardized coefficients"   
        path <- predict.l2boost(x, type = "coef")$coef.stand.path
      }
      else {
        if (is.null(ylab)) ylab <- "coefficients"   
        path <- predict.l2boost(x, type = "coef")$coef.path
      }
    }
    else {
      stop("type must be set to 'rho' or 'coef'\n")
    }
    # plot it
    if(is.null(xlim))xlim <- range(xval)
    if(is.null(ylim))ylim <- range(unlist(path))
    if(is.null(active.set)) active.set <- unique(l.crit)
    plot(range(xval), range(unlist(path)), type = "n",
         xlab = xlab,
         xlim = xlim,
         ylab = ylab,
         ylim = ylim, ...=...)
    #plot inactive set first
    if (type == "rho") {
      plot.lines(xval, setdiff(1:p, active.set), path, l.crit, FALSE, col)
      plot.lines(xval, active.set, path, l.crit, TRUE, col)
    }
    else {
      plot.lines(xval, active.set, path, l.crit, FALSE)
    }
  }
}
