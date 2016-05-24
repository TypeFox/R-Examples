#' Plot the prior, likelihood, and posterior on the same plot.
#' 
#' This function takes any object of class \code{Bolstad} and plots the prior,
#' likelihood and posterior on the same plot. The aim is to show the influence
#' of the prior, and the likelihood on the posterior.
#' 
#' 
#' @param x an object of class \code{Bolstad}.
#' @param \dots any other arguments to be passed to the \code{plot} function.
#' @note Note that \code{xlab}, \code{ylab}, \code{main}, \code{axes},
#' \code{xlim}, \code{ylim} and \code{type} are all used in the function so
#' specifying them is unlikely to have any effect.
#' @author James Curran
#' @keywords plots
#' @examples
#' 
#' # an example with a binomial sampling situation
#' results = binobp(4, 12, 3, 3, plot = FALSE)
#' decomp(results)
#' 
#' # an example with normal data
#' y = c(2.99,5.56,2.83,3.47)
#' results = normnp(y, 3, 2, 1, plot = FALSE)
#' decomp(results)
#' 
#' @export decomp
decomp = function(x, ...){
  if(class(x) != "Bolstad")
    stop("This function only works for objects of class Bolstad")
  
  oPar = par(mfrow = c(3, 1), mar = c(1, 1, 1, 1))
  with(x, {
    yLims = c(0, 1.1 * max(results$posterior, results$prior));
  
    plot(prior ~ param.x, ylim = yLims, type = "l", xlim = range(param.x),
         xlab = "", ylab = "", main = "",
         axes = FALSE, ...);
    polygon(param.x, prior, col = "red");
    box();
    r = legend("topleft", legend = "Prior",lty = 1, bty = "n", plot = FALSE)$text;
    text(r$x, r$y, "Prior", adj = 0);
    
    plot(likelihood ~ param.x, type="l",
         xlab = "", ylab = "", main = "",
         axes = FALSE, ...);
    polygon(param.x, likelihood, col = "green");
    box();
    r = legend("topleft", legend = "Prior",lty = 1, bty = "n", plot = FALSE)$text;
    text(r$x, r$y, "Likelihood", adj = 0);
    
    plot(posterior ~ param.x, ylim = yLims, type = "l",
         xlab = "", ylab = "", main = "",
         axes = F, ...);
    polygon(param.x, posterior, col = "blue");
    box();
    r = legend("topleft", legend = "Prior",lty = 1, bty = "n", plot = FALSE)$text;
    text(r$x, r$y, "Posterior", adj = 0);
  })
  
  par(oPar)
}
