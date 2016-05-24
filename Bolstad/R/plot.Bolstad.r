#' Plot method for objects of type Bolstad
#' 
#' A unified plotting method for plotting the prior, likelihood and posterior 
#' from any of the analyses in the book
#' 
#' The function provides a unified way of plotting the prior, likelihood and 
#' posterior from any of the functions in the library that return these 
#' quantities. It will produce an overlay of the lines by default, or separate 
#' panels if \code{overlay = FALSE}.
#' 
#' @param x A S3 object of class Bolstad
#' @param overlay if \code{FALSE} then up to three plots will be drawn 
#'   side-by-side
#' @param which Control which of the prior = 1, likelihood = 2, and posterior = 
#'   3, are plots. This is set to prior and posterior by default to retain 
#'   compatibility with the book
#' @param densCols The colors of the lines for each of the prior, likelihood and
#'   posterior
#' @param legendLoc The location of the legend, usually either \code{"topright"}
#'   or \code{"topleft"}
#' @param scaleLike If \code{TRUE}, then the likelihood will be scaled to have
#'   approximately the same maximum value as the posterior
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param main Title of plot
#' @param ylim Vector giving y coordinate range
#' @param cex Character expansion multiplier
#' @param \dots Any remaining arguments are fed to the \code{plot} command
#' @author James Curran
#' @keywords plot
#' @examples
#' 
#' x = rnorm(20,-0.5,1)
#' ## find the posterior density with a N(0,1) prior on mu
#' b = normnp(x,sigma=1)
#' plot(b)
#' plot(b, which = 1:3)
#' plot(b, overlay = FALSE, which = 1:3)
#' @export
plot.Bolstad = function(x, overlay = TRUE, which = c(1, 3), 
                        densCols = c("red","green","blue")[which],
                        legendLoc = "topleft", 
                        scaleLike = FALSE, 
                        xlab = eval(expression(x$name)), 
                        ylab = "", 
                        main = "Shape of prior and posterior", 
                        ylim = c(0, max(cbind(x$prior, x$likelihood, x$posterior)[,which]) * 1.1),
                        cex = 0.7, ...){
  
  which = sort(which)
  
  if(is.null(which) || length(which) <= 0 || length(which) > 3 || any(!grepl('^[1-3]+$', which))){
    stop("parameter which can only take vectors of length 3 containing the values 1, 2 and 3")
  }
  
  if(scaleLike){
    sf = max(x$posterior) / max(x$likelihood)
    x$likelihood = x$likelihood * sf
  }
  
  bLegend = !grepl("none", tolower(legendLoc))
  
  if(overlay){
    with(x,{
      Y = as.matrix(cbind(prior, likelihood, posterior)[,which]);
      
      plot(param.x, Y[,1], ylim = ylim, type="l",
           lty = (3:1)[which[1]], col = densCols[1],
           xlab = xlab, ylab = "",
           main = main, ...);
      
      i = 2;
      while(i <= ncol(Y)){
        lines(param.x, Y[,i], lty = (3:1)[which[i]], col = densCols[i]);
        i = i + 1;
      };
      
      if(bLegend){
        legend(legendLoc, lty = (3:1)[which], col = densCols, 
               legend = c("Prior", "Likelihood", "Posterior")[which], 
               bty = 'n', cex = cex);
      };
    })
  }else{
    oldpar = par(mfrow = c(1, length(which)), mai = c(0.7, 0.1, 0.2, 0.1), yaxs = 'i', xaxs = 'i')
    
    with(x,{
      Y = cbind(prior, likelihood, posterior)[,which]
    
      legend = c("Prior", "Likelihood", "Posterior")[which]
      
      plot(param.x, Y[,1], ylim = ylim, type="l",
           col = densCols[1],
           xlab = eval(expression(name)), ylab = "",
           main = legend[1], axes = FALSE, ...)
      axis(1)
      box()
      
      for(i in 2:ncol(Y)){
        plot(param.x, Y[,i], ylim = ylim, col = densCols[i], type = 'l', xlab = "", 
             main = legend[i], axes = FALSE, ...)
        box()
      }
      
      par(oldpar)
    })
  }
}
