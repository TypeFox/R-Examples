#' @title Plots to validate a Pareto Positive Stable (PPS) fit
#' @description Plots to validate a PPS fit (typically from \code{PPS.fit()}) with different comparisons between empirical and theoretical functions.
#' @param x a \code{PPSfit} Object.
#' @param which values from 1 to 4 indicating the type of plot.
#' @param ask an argument to control the plot window.
#' @param ylim optional argument to control the y limits of the histogram. It is included to prevent non-desired scales on the y-axis.
#' @param breaks optional argument to control the breakpoints of the histogram. See \code{hist} help for the details. It is included to prevent non-desired scales on the y-axis.
#' @param \dots other arguments.
#' 
#' @details The plots return:
#'     
#' 1. The histogram of the observations and the fitted PPS density (\code{which = 1}). Optional \code{ylim} and \code{breaks} arguments are provided to prevent frequent imbalances between density and histogram scales in real data: they work as the analogue arguments of the default \code{hist} function.
#' 2. The empirical distribution function of data and the cumulative distribution function of the fitted model (\code{which = 2}).
#' 3. A rank-size plot in log-log scale to check the Pareto or power-law behaviour of data (\code{which = 3}). In the X-axis the log of the observations appears; in the Y-axis, the log of the empirical survival function. If the scatter-plot is around a straight line, then the observations exhibit a power law behaviour. The plot also includes the curve with the theoretical survival function of the model specified in the first argument class \code{PPSfit}: only when \code{nu} is 1, that curve is going to be a straight line. 
#' 4. A plot in a double log-log scale to check the adequacy of data to the PPS model (\code{which = 4}). On one hand, the X-axis shows the double log of the observations divided by the scale parameter, while the Y-axis shows the log of minus the log of the empirical survival function. On the other hand, the straight line determined by the linear relation between the survival function and the scaled data in a double log-log scale, in relation to the argument class \code{PPSfit} is added. The proximity of the points in the scatter-plot to that straight line is an evidence in favour of a PPS behaviour of data.
#' 
#' @references Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
#' 
#' @seealso \code{\link{PPS.fit}}
#' 
#' @examples
#' data(forbes400)
#' fit <- PPS.fit(forbes400$NetWorth)
#' par(mfrow=c(2,2))
#' plot(fit)
#' dev.off()
#' plot(fit, which = 1, breaks = seq(0, 60, length.out = 60))

#' @importFrom graphics curve hist hist.default lines par plot
#' @importFrom grDevices dev.interactive
#' @importFrom stats ecdf

#' @export
plot.PPSfit <-
  function (x, which = 1:4, ask = prod(par("mfcol")) < length(which) && dev.interactive(), ylim, breaks, ...)
  {
    if (!class(x) == "PPSfit") 
      stop("Object must belong to class PPSfit")
    if (ask) {
      op <- par(ask = TRUE)
      on.exit(par(op))
    }
    par(mar = c(6, 4, 4, 2) + 0.1)
    show <- rep(FALSE, 4)
    show[which] <- TRUE
    lam <- x$estimate[[1]]
    if (is.null(x$Pareto)) x$Pareto <- FALSE
    if (x$Pareto==FALSE & !is.null(x$sigma)){
      sc <- x$sigma
      v <- x$estimate[[2]]
    }
    if (x$Pareto==FALSE & is.null(x$sigma)){   
      sc <- x$estimate[[2]]
      v <- x$estimate[[3]] 
    }  
    if (x$Pareto==TRUE & !is.null(x$sigma)){
      sc <- x$sigma
      v <- 1
    }
    if (x$Pareto==TRUE & is.null(x$sigma)){   
      sc <- x$estimate[[2]]
      v <- 1
    }
    obs <- x$obs
    obsName <- x$obsName
    if (missing(breaks)) breaks <- hist(obs, plot = FALSE, right = FALSE)$breaks
    else breaks <- hist(obs, plot = FALSE, right = FALSE, breaks)$breaks
    PPSDens <- function(y) dPPS(y, lam, sc, v)
    PPSProbs <- function(y) pPPS(sort(y), lam, sc, v)
    if (missing(ylim)) {
      ymax <- 1.06 * max(PPSDens(obs), na.rm = TRUE)
      ylim <- c(0, ymax)
    }
    if (show[1]) {
      hist.default(obs, right = FALSE, freq = FALSE, 
                   ylim = ylim, 
                   main=NULL,
                   xlab="x",
                   breaks, 
                   ...)       
      curve(PPSDens, min(breaks) - 1, max(breaks) + 1, add = TRUE, 
            ylab = NULL,col=2,lty=2,lwd=2)
      #legend("topright", legend = substitute(paste("PPS(",lambda, " = ", nn1 ,", ", sigma, " = ", nn2, ", ", nu, " = ", nn3,")"), list(nn1=round(lam,3),nn2=round(sc,3),nn3=round(v,3))), bty = "n", col = 2, lty = 2, lwd = 2)
    }
    if (show[2]) {
      ecdf2<-function(y) length(y)*ecdf(y)(y)/(length(y)+1)
      plot(obs, ecdf2(obs), 
           xlab="x",
           ylab="Cumulative probability", 
           ...)
      lines(sort(obs), PPSProbs(sort(obs)),col=2,lty=2,lwd=2)
      #legend("bottomright",legend=substitute(paste("PPS(",lambda, " = ", nn1 ,", ", sigma, " = ", nn2, ", ", nu, " = ", nn3,")"), list(nn1=round(lam,3),nn2=round(sc,3),nn3=round(v,3))), bty = "n", col = 2, lty = 2, lwd = 2)
    }
    if (show[3]) {
      rango<-function(y) -(length(y)*ecdf(y)(y)-(length(y)+1))
      plot(log(obs), log(rango(obs)), 
           xlab="log(x)",
           ylab="log(rank(x))",
           ...)
      lines(log(sort(obs)), log(x$n*(1-PPSProbs(sort(obs)))),col=2,lty=2,lwd=2)
      #legend("bottomleft", legend = substitute(paste("PPS(",lambda, " = ", nn1 ,", ", sigma, " = ", nn2, ", ", nu, " = ", nn3,")"), list(nn1=round(lam,3),nn2=round(sc,3),nn3=round(v,3))), bty = "n", col = 2, lty = 2, lwd = 2)
    }
    if (show[4]) {
      rango<-function(y) -(length(y)*ecdf(y)(y)-(length(y)+1))
      if (min(obs) < sc+0.000001) obs <- obs[which(obs >= sc+0.000001)]
      plot(log(log(obs/sc)), log(-log(rango(obs)/(x$n+1))), 
           xlab="log(log(x/scale))",
           ylab="log(-log(rank(x)/(n+1)))",
           ...)
      lines(log(log(sort(obs)/sc)), log(-log(1-PPSProbs(sort(obs)))),col=2,lty=2,lwd=2)
      #legend("topleft",legend=substitute(paste("PPS(",lambda, " = ", nn1 ,", ", sigma, " = ", nn2, ", ", nu, " = ", nn3,")"), list(nn1=round(lam,3),nn2=round(sc,3),nn3=round(v,3))), bty = "n", col = 2, lty = 2, lwd = 2)
    }
    invisible()
  }