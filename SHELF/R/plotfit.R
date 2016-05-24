#' Plot the fitted density function for one or more experts
#' 
#' Plots the fitted density function for one or more experts. Can also plot a
#' fitted linear pool if more than one expert. If plotting the density function
#' of one expert, or the linear pool only, can also indicated desired lower and
#' upper fitted quantiles.
#' 
#' 
#' @param fit The output of a \code{fitdist} command.
#' @param d The distribution fitted to each expert's probabilities. Options are
#' \code{"normal"}, \code{"t"}, \code{"gamma"}, \code{"lognormal"},
#' \code{"logt"},\code{"beta"}, \code{"hist"} (for a histogram fit), and
#' \code{"best"} (for best fitting)
#' @param int Set \code{int = TRUE} to use interactive plotting (using the
#' shiny package). If plotting for a single expert, the argument \code{d} is
#' ignored, as distributions can be chosen within the display. If plotting for
#' multiple experts, feedback quantiles are not displayed, and the argument
#' \code{lp} is ignored, as the option to show a linear pool can be chosen
#' within the display.
#' @param xl The lower limit for the x-axis. The default is the 0.001 quantile
#' of the fitted distribution (or the 0.001 quantile of a fitted normal
#' distribution, if a histogram fit is chosen).
#' @param xu The upper limit for the x-axis. The default is the 0.999 quantile
#' of the fitted distribution (or the 0.999 quantile of a fitted normal
#' distribution, if a histogram fit is chosen).
#' @param ql A lower quantile to be indicated on the density function plot.
#' Only displayed when plotting the density function for a single expert.
#' @param qu An upper quantile to be indicated on the density function plot.
#' Only displayed when plotting the density function for a single expert.
#' @param lp For multiple experts, set \code{lp=TRUE} to plot a linear pool.
#' @param ex If judgements have been elicited from multiple experts, but a
#' density plot for one expert only is required, the expert to be used in the
#' plot.
#' @param sf The number of significant figures to be displayed for the
#' parameter values.
#' @param ind If plotting a linear pool, set \code{ind=FALSE} to suppress
#' plotting of the individual density functions.
#' @param lpw A vector of weights to be used in linear pool, if unequal
#' weighting is desired.
#' @author Jeremy Oakley <j.oakley@@sheffield.ac.uk>
#' @examples
#' 
#' \dontrun{
#' # Two experts
#' # Expert 1 states P(X<30)=0.25, P(X<40)=0.5, P(X<50)=0.75
#' # Expert 2 states P(X<20)=0.25, P(X<25)=0.5, P(X<35)=0.75
#' # Both experts state 0<X<100. 
#' 
#' v <- matrix(c(30, 40, 50, 20, 25, 35), 3, 2)
#' p <- c(0.25, 0.5, 0.75)
#' myfit <- fitdist(vals = v, probs = p, lower = 0, upper = 100)
#' 
#' # Plot both fitted densities, using the best fitted distribution
#' plotfit(myfit)
#' 
#' # Plot a fitted beta distribution for expert 2, and show 5th and 95th percentiles
#' plotfit(myfit, d = "beta", ql = 0.05, qu = 0.95, ex = 2)
#' 
#' # Use interactive plotting for for expert 2, and show 5th and 95th percentiles
#' plotfit(myfit, int = T, ex = 2)
#' 
#' # Plot a linear pool, giving double weight to expert 1
#' plotfit(myfit,  lp = T, lpw = c(2,1))
#' 
#' # Use interactive plotting, giving double weight to expert 1, if a linear pool is displayed
#' plotfit(myfit,  int = T, lpw = c(2,1))
#' 
#' # Plot a linear pool, giving double weight to expert 1, 
#' # show 5th and 95th percentiles, supress plotting of individual distributions, 
#' # and force use of Beta distributions
#' plotfit(myfit, d = "beta",  lp = T, lpw = c(2,1), ql = 0.05, qu = 0.95, ind=FALSE )
#' }
#' @import graphics
#' @export
plotfit <-
function(fit, d = "best", int = FALSE, xl = -Inf, xu = Inf, ql = NA, qu = NA, lp = FALSE, ex = NA, sf = 3, ind = TRUE, lpw = 1){

  if(d=="beta" & (min(fit$limits) == -Inf | max(fit$limits) == Inf )){stop("Parameter limits must be finite to fit a beta distribution")}
  if(d=="gamma" & min(fit$limits) < 0 ){stop("Lower parameter limit must be non-negative to fit a gamma distribution")}
  if(d=="lognormal" & min(fit$limits) < 0 ){stop("Lower parameter limit must be non-negative to fit a log normal distribution")}
  if(d=="logt" & min(fit$limits) < 0 ){stop("Lower parameter limit must be non-negative to fit a log t distribution")}
  if(is.na(ql)==F & (ql <0 | ql>1 )){stop("Lower feedback quantile must be between 0 and 1")}
  if(is.na(qu)==F & (qu <0 | qu>1 )){stop("Upper feedback quantile must be between 0 and 1")}
  
  
  if(nrow(fit$vals)>1 & is.na(ex)==T & lp==F){
    if(xl == -Inf & min(fit$limits[,1]) > -Inf){xl <- min(fit$limits[,1]) }
    if(xu == Inf & max(fit$limits[,2]) < Inf){xu <- max(fit$limits[,2]) }
    if(int == FALSE){plotgroup(fit, xl, xu, d, bw = T)}else{
      shinyplotgroup(fit, xl, xu, lpw)
    }
  }
  
  if(nrow(fit$vals)>1 & lp==T){
    if(xl == -Inf & min(fit$limits[,1]) > -Inf){xl <- min(fit$limits[,1]) }
    if(xl == -Inf & min(fit$limits[,1]) == -Inf){
      f1 <- feedback(fit, quantiles=0.01, dist=d)
      xl <- min(f1$expert.quantiles)
    }
    
    if(xu == Inf & max(fit$limits[,2]) < Inf){xu <- max(fit$limits[,2]) }
    
    if(xu == Inf & max(fit$limits[,2]) == Inf){
      f2 <- feedback(fit, quantiles=0.99, dist=d)
      xu <- max(f2$expert.quantiles)
    }
    if(int == FALSE){plotlinearpool(fit, xl, xu, ql, qu , d, ind, lpw)}else{
      shinyplotgroup(fit, xl, xu, lpw)
    }
    
  }
  
  if(nrow(fit$vals)>1 & is.na(ex)==F){
    if(xl == -Inf & fit$limits[ex,1] > -Inf){xl <- fit$limits[ex,1] }
    if(xu == Inf & fit$limits[ex,2] < Inf){xu <- fit$limits[ex,2] }
    if(int == FALSE){plotsingle(fit, d, xl, xu, ql, qu, sf, ex)}else{
      shinyplotsingle(fit, xl, xu, ql, qu, ex)
    }
    
  }
  
  if(nrow(fit$vals)==1){
    if(int == FALSE){plotsingle(fit, d, xl, xu, ql, qu, sf, ex = 1)}else{
      shinyplotsingle(fit, xl, xu, ql, qu, ex = 1)
    }
    
  }	
}
