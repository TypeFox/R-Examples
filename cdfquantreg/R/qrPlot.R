#' @title Plot Fitted Values/Residuals of A Cdfqr Object or Distribution
#' @aliases plot.cdfqr
#' @description Plot Fitted Values/Residuals of A Cdfqr Object or Distribution
#' @param x If the plot is based on the fitted values, provide a fitted cdfqr object.
#'
#' @param mu,sigma, fd, sd alternatively, mu and sigma, and the distribution can be specified
#' @param n The number of random variates to be generated for user specified plot. 
#' @param type Currently only fitted values are available for generating plots.
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#' @param ... other plot parameters pass onto \code{\link[graphics]{plot}}.
#' @examples
#' data(cdfqrExampleData)
#' fit <- cdfquantreg(crc99 ~ vert | confl, 't2','t2', data = JurorData)
#' plot(fit)
#' 
#' 
#' @method plot cdfqr
#' @export
#' @import graphics
#' @importFrom MASS truehist

plot.cdfqr <- function(x, mu = NULL, sigma = NULL, fd = NULL, sd = NULL, n = 10000, 
                       type = c("fitted"), ...) {
  
  # If plot based on the fitted model
  if(is.null(mu)){
    ydata <- x$y
    fd <- x$family$fd
    sd <- x$family$sd
    mu <- mean(x$fitted$mu)
    sigma <- exp(mean(log(x$fitted$sigma)))
    # smooth the boundary values
    fitted <- scaleTR(x$fitted$full)
    fit_d <- dq(fitted, mu, sigma, fd, sd)
    xtem <- data.frame(x = fitted, y = as.numeric(fit_d))
    xtem <- xtem[order(xtem$x),]
    
    par(mfrow=c(1,2))
    # plot the density
#     if(missing(xlim)) xlim <- c(0,1)
#     if(missing(ylim)) ylim <- c(0, max(max(c(xtem$x, xtem$y))))
    MASS::truehist(ydata, col = "white", ymax = max(xtem$y) + 0.1,
                   main = "Data (histogram) \n fitted by model (line)",
                   ...)
    graphics::lines(xtem$x, xtem$y, lty = 1, lwd = 2)

    
    # plot data against fitted values
    plot(ydata,x$fitted$full, xlab = 'Observations', ylab = 'Fitted',
         main = "Observations vs. Fitted",type ='p', pch=20)
    graphics::abline(0, 1)
    
  }else{
    #based on input mu and sigma, and given distribution
   ydata <- as.numeric(rq(n, mu, sigma, fd, sd))
   fit_d <- as.numeric(dq(ydata, mu, sigma, fd, sd))
   xtem <- data.frame(x = ydata, y =fit_d)
   xtem <- xtem[order(xtem$x),]
   xlim <- c(0,1)
   # ylim <- c(0, max(max(c(xtem$x, xtem$y))))
   MASS::truehist(ydata, col = "white", ymax = max(xtem$y) + 0.1, ...)
   graphics::lines(xtem$x, xtem$y, lty = 1, lwd = 2)
  }
  
  invisible()
}

