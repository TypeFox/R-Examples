#' @title Plot gravity model
#' @description Diagnostic plots gravity model with 6 optional plots.
#' 
#' @param  x         Object of class gravity 
#' @param  type      Type of plot (default 1, model structure I)
#' @param ...        Ignored
#'
#' @return defined plot 
#'
#' @note Plot types available: 1 - Model structure I, 2 - Model structure II, 3 - Q-Q Normal - Origin random effects, 4 - Q-Q Normal - Residuals , 5 - Fitted values, 6 - Distribution of observed verses predicted		
#' @note Depends: nlme, lattice
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org> and Melanie Murphy <melanie.murphy@@uwyo.edu>
#'
#' @references
#' Murphy, M. A. & J.S. Evans. (in prep). "GenNetIt: gravity analysis in R for landscape genetics" 
#' @references
#' Murphy M.A., R. Dezzani, D.S. Pilliod & A.S. Storfer (2010) Landscape genetics of high mountain frog metapopulations. Molecular Ecology 19(17):3634-3649 
#'  
#' @method plot gravity 
#' @export
plot.gravity <- function(x, type = 1, ...) {
  if(type == 1) {
  # MODEL STRUCTURE I  
  	graphics::plot(stats::fitted(x$gravity, level=0), x$y, xlab = "Fitted Values (DPS)",
                   ylab="Observed Values", main="Model Structure (I)", pch=16, ...)
                   graphics::abline(0, 1, col = "blue")
    }
  if(type == 2) {
  # MODEL STRUCTURE II
  	stats::scatter.smooth(stats::fitted(x$gravity), stats::residuals(x$gravity, type="pearson"), 
                    ylab="Innermost Residuals", main="Model Structure (II)",
					xlab="Fitted Values", pch=16)
                    graphics::abline(h = 0, col = "red")
	}
  if(type == 3) {			  
  # Q-Q NORMAL - ORIGIN RANDOM EFFECTS		   
      stats::qqnorm(nlme::ranef(x$gravity)[[1]], main="Q-Q Normal - Origin Random Effects", pch=16)
             stats::qqline(nlme::ranef(x$gravity)[[1]], col="red")
    }
  if(type == 4) {
  # Q-Q NORMAL - RESIDUALS		   
  	stats::qqnorm(stats::residuals(x$gravity, type="pearson"), main="Q-Q Normal - Residuals", pch=16)
                  stats::qqline(stats::residuals(x$gravity, type="pearson"), col="red")
    }  
  if(type == 5) {
  # FITTED VALUES	
    options(warn=-1)	   
  	graphics::boxplot(stats::residuals(x$gravity, type="pearson", level=1) ~ x$groups,
              ylab="Innermost Residuals", xlab="Origin",
              notch=T, varwidth = T, at=rank(nlme::ranef(x$gravity)[[1]]))
              graphics::axis(3, labels=format(nlme::ranef(x$gravity)[[1]], dig=2), cex.axis=0.8,
  			  at=rank(nlme::ranef(x$gravity)[[1]]))
              graphics::abline(h=0, col="darkgreen")  
    }  
  if(type == 6) {
  # DISTRIBUTION OF OBSERVED VS. PRED		
    oden <- stats::density(x$y)
    pden <- stats::density(stats::predict(x$gravity)) 
    graphics::plot(oden, type="n", main="", xlim=c(min(x$y), max(x$y)), 
         ylim=c(min(oden$y,pden$y), max(oden$y,pden$y)))     
            graphics::polygon(oden, col=grDevices::rgb(1,0,0,0.5))
            graphics::polygon(pden, col=grDevices::rgb(0,0,1,0.5))
      graphics::legend("topright", legend=c("Obs","Pred"), 
	              fill=c(grDevices::rgb(1,0,0,0.4), grDevices::rgb(0,0,1,0.4)))
    }
}
