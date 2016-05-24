#' @title Summarizing Bayesian change point analysis results
#' @description
#' Summary and print methods for class \code{bcp}.
#'
#' @param x the result of a call to \code{bcp()}.
#' @param object the result of a call to \code{bcp()}.
#' @param digits the number of digits displayed in the summary statistics.
#' @param ... (optional) additional arguments, ignored.
#' 
#' @author Xiaofei Wang, Chandra Erdman, and John W. Emerson
#' @return 
#' The matrix of results is returned invisibly.
#' @details 
#' The functions print (and return invisibly) the estimated posterior probability of a change point for each position and the estimated posterior means. These results are modeled after the summary method of the \code{coda} package (Plummer \emph{et al.}, 2006). If \code{return.mcmc=TRUE} (i.e., if full MCMC results are returned), \code{bcp} objects can be converted into \code{mcmc} objects to view \code{mcmc} summaries -- see examples below.
#' @seealso \code{\link{bcp}} and \code{\link{plot.bcp}}.
#' @examples 
#' ##### A random sample from a few normal distributions #####
#' testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
#' bcp.0 <- bcp(testdata)
#' summary(bcp.0)
#' plot(bcp.0, main="Univariate Change Point Example")
#' 
#' ##### An MCMC summary from the ``coda'' package #####
#' if (require("coda")) {
#'   bcp.0 <- bcp(testdata, return.mcmc=TRUE)
#'   bcp.mcmc <- as.mcmc(t(bcp.0$mcmc.means))
#'   summary(bcp.mcmc)
#'   heidel.diag(bcp.mcmc) # an example convergence diagnostic
#'   # from the coda package.
#' }
#' 
#' @keywords datasets
#' @export
summary.bcp <-
  function (object, digits = max(3, .Options$digits - 3), ...) {
  
  if (is.null(colnames(object$data))) {
    colnames(object$data) <- c("Loc", paste(rep("X", ncol(object$data)-1), 
                                   1:(ncol(object$data)-1), sep=""))
  }

  statnames <- c("Probability", colnames(object$posterior.mean))
  varstats <- matrix(NA, nrow = max(object$data[,1]), 
                     ncol = length(statnames), 
                     dimnames = list(1:max(object$data[,1]), statnames))
						   
  varstats[,1] <- object$posterior.prob
  varstats[,2:ncol(varstats)] <- object$posterior.mean

  cat("\nBayesian Change Point (bcp) summary:\n\n")
  cat("\nProbability of a change in mean and posterior means:\n\n")
  print(varstats, digits=digits)
  cat("\n")

  invisible(varstats)
}

#' @rdname summary.bcp
#' @export
print.bcp <- function(x, digits = max(3, .Options$digits - 3), ...) {
  summary(x, digits=digits, ...)
}
#' @title Estimate the probability of a change point in a specified interval
#' 
#' @description
#' The function \code{interval.prob()} estimates the probability of at least one
#' change point in the specified interval of sequential observations; it may only be used when \code{return.mcmc=TRUE}.
#'
#' @param object the result of a call to \code{bcp()}.
#' @param start the starting index of the interval.
#' @param end the ending index of the interval.
#' 
#' @author Xiaofei Wang, Chandra Erdman, and John W. Emerson
#' @details For sequential data only, the function returns an estimate of the posterior probability of at least one change point in the specified interval.
#' @note \code{return.mcmc} must be \code{TRUE}.
#' @seealso \code{\link{bcp}} and \code{\link{plot.bcp}}.
#' @examples 
#' ##### A random sample from a few normal distributions #####
#' testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
#' bcp.0 <- bcp(testdata, return.mcmc=TRUE)
#' plot(bcp.0, main="Univariate Change Point Example")
#' interval.prob(bcp.0, 45, 55)
#' 
#' @keywords datasets
#' @export	
interval.prob <- function(object, start, end) {
  if (attr(object, "structure") == "graph") stop("Method is not implemented for graph change points.")
  if (!object$return.mcmc) stop("bcp must be run with return.mcmc=TRUE")
  return( sum(apply(object$mcmc.rhos[start:end,-c(1:object$burnin)], 2, sum) > 0) /
          object$mcmc)
}

#' @title Plotting univariate Bayesian change point results
#' @aliases legacyplot
#' @description
#' \code{legacyplot()} produces summary plots of the results of \code{bcp()} when used for univariate analysis; it was the default method prior to package version 3.0.0.
#' @param x the result of a call to \code{bcp()}.
#' @param ... (optional) additional arguments, ignored.
#' 
#' @author Chandra Erdman and John W. Emerson
#' @details \code{legacyplot()} produces the following plots using \code{base} graphics:
#'
#' Posterior Means: location in the sequence versus the posterior mean over the iterations. 
#' 
#' Posterior Probability of a Change: location in the sequence versus the relative frequency of iterations which resulted in a change point. 
#' 
#' @seealso \code{\link{plot.bcp}}, \code{\link{bcp}}, \code{\link{summary.bcp}}, and \code{\link{print.bcp}} for complete results and summary statistics.
#' @examples 
#' ##### A random sample from a few normal distributions #####
#' testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
#' bcp.0 <- bcp(testdata, return.mcmc=TRUE)
#' legacyplot(bcp.0)
#' 
#' @keywords datasets
#' @export
legacyplot <- function(x, ...) { 	
  if ((is.matrix(x$data) && ncol(x$data)>2 ) && 
       !(attr(x, "model") == "multivariate" || attr(x, "structure") == "series"))
    stop("Legacy bcp plot invalid for multivariate bcp object or graph bcp object.")
  posterior.prob <- x$posterior.prob
  posterior.prob[length(posterior.prob)] <- 0
		
  op <- par(mfrow=c(2,1),col.lab="black",col.main="black")
  op2 <- par(mar=c(0,4,4,2),xaxt="n", cex.axis=0.75)
  plot(1:nrow(x$data), x$data[,2], col="grey", 
       pch=20, xlab="", ylab="Posterior Mean", 
       main="Posterior Means and Probabilities of a Change", ...)
  lines(x$posterior.mean, lwd=2)
  par(op2)
  op3 <- par(mar=c(5,4,0,2), xaxt="s", cex.axis=0.75)
  plot(1:length(x$posterior.mean), posterior.prob, 
       yaxt="n", type="l", ylim=c(0,1),
       xlab="Location", ylab="Posterior Probability", main="")
  axis(2, yaxp=c(0, 0.9, 3))
  par(op3)
  par(op)
}

#' @title Plotting Bayesian change point results
#' 
#' @description
#' \code{plot.bcp()} produces summary plots of the results of \code{bcp()}. Currently, only the summary plots for serial data are implemented. 
#' If an adjacency structure (adj) is provided, then the data are assumed to reside on nodes of a general graph. Additional parameters are used in this graph change point model.
#'
#' @details 
#' \code{plot.bcp()} produces the following plots using \code{grid} graphics instead of \code{base}:
#'   
#'  Posterior Means: location in the sequence versus the posterior means over the iterations. 
#' 
#'  Posterior Probability of a Change: location in the sequence versus the relative frequency of iterations which resulted in a change point. 
#' 
#' @param x the result of a call to \code{bcp()}.
#' @param separated logical. If set to \code{TRUE} and the data is multivariate, each series is plotted separately.
#' @param outer.margins (optional) list of units specifying the left, bottom, right and top margins.  For more information on units, see the documentation for \code{grid}.
#' @param lower.area (optional) unit specifying the proportion of the plot occupied by the posterior probabilities of change points.
#' @param size.points (optional) unit specifying the size of the data points.
#' @param pch.points (optional) unit specifying the style of the data points.
#' @param colors (optional) vector specifying the colors in which to plot each data series.
#' @param main (optional) plot title. Use \code{""} for no title.
#' @param xlab (optional) a character string specifying the x-axis label. Defaults to "Location".
#' @param xaxlab (optional) a vector having length equal to the number of observations giving the x-axis tick labels. Defaults to the sequence from 1 to n.
#' @param cex.axes (optional) list specifying the sizes of the axes labels. \code{cex.xaxis} specifies the size of the x-axis label, \code{cex.yaxis.lower} specifies the size of the y-axis label of the posterior probability plot, \code{cex.yaxis.upper.default} specifies the size of the y-axis labels of the posterior means plot when the series are displayed in a single plot, and \code{cex.yaxis.upper.separated} specifies the size of the y-axis labels of the posterior means plots when each series is plotted separately.
#' @param ... (optional) additional arguments, ignored.
#' 
#' @author Xiaofei Wang, Chandra Erdman, and John W. Emerson
#' 
#' @seealso \code{\link{legacyplot}}, \code{\link{bcp}}, \code{\link{summary.bcp}}, and \code{\link{print.bcp}} for complete results and summary statistics.
#' @examples 
#' testdata <- cbind( c(rnorm(50), rnorm(50, -5, 1), rnorm(50)),
#' c(rnorm(50), rnorm(50, 10.8, 1), rnorm(50, -3, 1)) )
#' bcp.0 <- bcp(testdata)
#' plot(bcp.0, main="Multivariate (k=2) Change Point Example")
#' plot(bcp.0, separated=TRUE, main="Multivariate (k=2) Change Point Example")
#' 
#' @keywords datasets
#' @export
plot.bcp <- function(x, separated = FALSE, 
                     outer.margins = list(left=unit(4, "lines"),
                                          bottom=unit(3, "lines"),
                                          right=unit(2, "lines"), 
                                          top=unit(2, "lines")),
                     lower.area = unit(0.33, "npc"),
                     size.points = unit(0.25, "char"),
                     pch.points = 20,
                     colors = NULL,
                     main = NULL,
                     xlab = NULL,
                     xaxlab = NULL,
                     cex.axes = list(cex.xaxis = 0.75,
                                     cex.yaxis.lower = 0.75,
                                     cex.yaxis.upper.default = 0.75,
                                     cex.yaxis.upper.separated = 0.5),
                     ...) {
  # if (!require(grid)) stop("library(grid) is required and unavailable.\n\n")
  if (attr(x, "structure") == "graph")   
    stop("Plotting functionality for data on a graph has not yet been implemented.")

  grid.newpage()

  thisjust <- c("left", "bottom")
  n <- length(x$posterior.prob)
  if (!is.matrix(x$data)) {
    x$data <- matrix(x$data, ncol=1) 
    x$posterior.mean <- matrix(x$posterior.mean, ncol=1)
  }
  m <- ncol(x$data)
  if (is.null(colors)) {
    colors <- 2:(m+1)
  } else {
    if (length(colors)!=m) stop("length(colors) must equal number of series.")
  }
  if (is.null(main)) main <- "Posterior Means and Probabilities of a Change"
  if (is.null(xlab)) xlab <- "Location"
  
  vp.main <- viewport(x=outer.margins$left, y=outer.margins$bottom,
                      width=unit(1, "npc")-outer.margins$right-outer.margins$left,
                      height=unit(1, "npc")-outer.margins$top-outer.margins$bottom,
                      just=thisjust, name="main", clip="off")
  pushViewport(vp.main)
  grid.rect()
  grid.text(main, 0.5, unit(1, "npc")+unit(1,"lines"), 
            gp=gpar(fontface="bold"))

  # Upper plot
  if (attr(x, "model") == "multivariate") {
    if (!separated) {
      pushViewport(viewport(x=unit(0, "npc"), y=lower.area,
                            width=unit(1, "npc"), height=unit(1, "npc") - lower.area,
                            just=thisjust, name="upper", clip="off",
                            default.units="native",
                            xscale=range(pretty(x$data[,1])),
                            yscale=range(pretty(x$data[,-1], n=10))))
      grid.rect()
      grid.yaxis(gp=gpar(cex=cex.axes$cex.yaxis.upper.default), main=FALSE)
      grid.text("Posterior Means", unit(-3, "lines"), 0.5, rot=90)
      for (i in 2:m) {
        grid.points(x$data[,1], x$data[,i], size=size.points, gp=gpar(col=colors[i]-1),
                    pch=pch.points)
        grid.lines(unique(x$data[,1]), x$posterior.mean[,i-1], gp=gpar(col=colors[i]-1),
                   default.units="native")
      }
      popViewport(1)
    } else {
      pushViewport(viewport(x=unit(0, "npc"), y=lower.area,
                            width=unit(1, "npc"), height=unit(1, "npc") - lower.area,
                            just=thisjust, name="upper", clip="off"))
      grid.text("Posterior Means", unit(-3, "lines"), 0.5, rot=90)
      yloc <- FALSE
      for (i in 2:m) {
        pushViewport(viewport(x=unit(0, "npc"),
                              y=unit((i-1)/m, "npc"),
                              width=unit(1, "npc"),
                              height=unit(1/m, "npc"),
                              just=thisjust, name="upper", clip="off",
                              default.units="native",
                              xscale=range(pretty(x$data[,1])),
                              yscale=range(pretty(x$data[,i]))))
        grid.rect()
        grid.yaxis(gp=gpar(cex=cex.axes$cex.yaxis.upper.separated), main=yloc)
        grid.points(x$data[,1], x$data[,i], size=size.points, gp=gpar(col=colors[i]-1),
                    pch=pch.points)
        grid.lines(unique(x$data[,1]), x$posterior.mean[,i-1], gp=gpar(col=colors[i]-1),
                      default.units="native")
        popViewport(1)
        yloc <- !yloc
      }
      popViewport(1)
    }
  } else {
    if (!separated) {
      pushViewport(viewport(x=unit(0, "npc"), y=lower.area,
                            width=unit(1, "npc"), height=unit(1, "npc") - lower.area,
                            just=thisjust, name="upper", clip="off",
                            default.units="native",
                            xscale=range(pretty(x$data[,1])),
                            yscale=range(pretty(x$data[,2]))))
      grid.rect()
      grid.yaxis(gp=gpar(cex=cex.axes$cex.yaxis.upper.default), main=FALSE)
      grid.text("Posterior Means", unit(-3, "lines"), 0.5, rot=90)
  #    for (i in 1:m) {
        grid.points(x$data[,1], x$data[,2], size=size.points, gp=gpar(col=colors[1]),
                    pch=pch.points)
        grid.lines(1:n, x$posterior.mean[,1], gp=gpar(col=colors[1]),
                   default.units="native")
  #    }
      popViewport(1)
    } else {
      pushViewport(viewport(x=unit(0, "npc"), y=lower.area,
                            width=unit(1, "npc"), height=unit(1, "npc") - lower.area,
                            just=thisjust, name="upper", clip="off"))
      grid.text("Posterior Means", unit(-3, "lines"), 0.5, rot=90)
      yloc <- TRUE
      m <- m - 1
      for (i in 2:(m+1)) {
        pushViewport(viewport(x=unit(0, "npc"),
                              y=unit((m-i+1)/m, "npc"),
                              width=unit(1, "npc"),
                              height=unit(1/m, "npc"),
                              just=thisjust, name="upper", clip="off",
                              default.units="native",
                              xscale=range(pretty(x$data[,1])),
                              yscale=range(pretty(x$posterior.mean[,i]))))
        grid.rect()
        grid.yaxis(gp=gpar(cex=cex.axes$cex.yaxis.upper.separated), main=yloc)
        if (yloc) {
          grid.text(colnames(x$posterior.mean)[i], x=unit(-2, "lines"), rot=90, gp=gpar(cex.xaxis=0.5))
        } else {
          grid.text(colnames(x$posterior.mean)[i], x=1.05, rot=90, gp=gpar(cex.xaxis=0.5))
        }

        grid.points(unique(x$data[,1]), x$posterior.mean[,i], size=size.points, gp=gpar(col=colors[1]),
                    pch=pch.points)
        grid.lines(unique(x$data[,1]), x$posterior.mean[,i], gp=gpar(col=colors[1]),
                      default.units="native")
        popViewport(1)
        yloc <- !yloc
      }
      popViewport(1)
    }
  }

  # Lower plot
  pushViewport(viewport(x=unit(0, "npc"), y=unit(0, "npc"),
                        width=unit(1, "npc"), height=lower.area,
                        just=thisjust, name="lower", clip="off",
                        default.units="native",
                        xscale=range(pretty(x$data[,1])),
                        yscale=c(-0.05, 1.05)))
  grid.rect()
  grid.yaxis(gp=gpar(cex=cex.axes$cex.yaxis.lower))
  if (!is.null(xaxlab)) {
    grid.xaxis(at = 1:n, label=xaxlab, gp=gpar(cex=cex.axes$cex.xaxis))
  } else grid.xaxis(gp=gpar(cex=cex.axes$cex.xaxis))
  grid.text("Posterior Probability", unit(-3, "lines"), 0.5, rot=90)
  grid.text(xlab, 0.5, unit(-2, "lines"))
  grid.lines(1:max(x$data[,1]), x$posterior.prob, default.units="native")
  popViewport(2) # vp.main and vp.lower

}

#' @title Extract model residuals
#' 
#' @description
#' residuals method for class \code{bcp}.
#'
#' @param object the result of a call to \code{bcp()}.
#' @param ... (optional) additional arguments, ignored.
#' 
#' @author Xiaofei Wang, Chandra Erdman, and John W. Emerson
#' @return 
#' Residuals extracted from the \code{bcp} object.
#' @seealso \code{\link{bcp}} and \code{\link{plot.bcp}}
#' @examples 
#' ##### A random sample from a few normal distributions #####
#' testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
#' bcp.0 <- bcp(testdata)
#' residuals(bcp.0)
#' 
#' @keywords datasets
#' @export
residuals.bcp <- function(object, ...) { 
  if (attr(object, "model") == "multivariate")
    return(object$data[,-1] - fitted.bcp(object))
  else 
    return(object$data[,2] - fitted.bcp(object))
}
#' @title Extract model fitted values
#' 
#' @description
#' fitted method for class \code{bcp}.
#'
#' @param object the result of a call to \code{bcp()}.
#' @param ... (optional) additional arguments, ignored.
#' 
#' @author Xiaofei Wang, Chandra Erdman and John W. Emerson
#' @return 
#' Fitted values extracted from the \code{bcp} object.
#' @seealso \code{\link{plot.bcp}}, \code{\link{summary.bcp}}, and \code{\link{print.bcp}} for summaries of the results.
#' @examples 
#' ##### A random sample from a few normal distributions #####
#' testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
#' bcp.0 <- bcp(testdata)
#' residuals(bcp.0)
#' 
#' @keywords datasets
#' @export
fitted.bcp <- function(object, ...) {
  if (attr(object, "model") == "multivariate")
    return(object$posterior.mean[object$data[,1],])
  else
    return(object$posterior.mean[object$data[,1],1])
}


