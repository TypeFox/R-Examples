#' @title Gravity model
#' @description Implements Murphy et al., (2010) gravity model
#' 
#' @param y              Name of dependent variable
#' @param x              Character vector of independent variables
#' @param d              Name of column contaning distance
#' @param group          Name of grouping column (from or to)
#' @param data           data.frame object containing model data
#' @param ln             Natural log transform data (TRUE/FALSE) 
#' @param constrained    Specify constrained model, if FALSE a linear model (lm) is run (TRUE/FALSE)
#' @param ...            Additional argument passed to nlme or lm      
#' 
#' @return formula         Model formula  
#' @return gravity         Gravity model
#' @return AIC             AIC value for selected model
#' @return x               data.frame of independent variables
#' @return y               Vector of dependent variable
#' @return groups          Ordered factor vector of grouping variable
#' @return fit             Model Fitted Values
#'
#' @note The "group" factor defines the singly constrained direction (from or to) and the grouping structure for the origins.
#'
#' @note Depends: nlme, lattice
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org> and Melanie Murphy <melanie.murphy@@uwyo.edu>
#'
#' @references
#' Murphy, M. A. & J.S. Evans. (in prep). "GenNetIt: gravity analysis in R for landscape genetics" 
#' @references
#' Murphy M.A., R. Dezzani, D.S. Pilliod & A.S. Storfer (2010) Landscape genetics of high mountain frog metapopulations. Molecular Ecology 19(17):3634-3649 
#'
#' @examples 
#' library(nlme)
#' data(ralu.model)
#' 
#' # Gravity	
#' x = c("DISTANCE", "DEPTH_F", "HLI_F", "CTI_F", "cti", "ffp")
#' ( gm <- gravity(y = "DPS", x = x, d = "DISTANCE", group = "FROM_SITE", 
#'                 data = ralu.model, ln = FALSE) )
#' 
#' # Plot gravity results
#'  par(mfrow=c(2,3))
#'    for (i in 1:6) { plot(gm, type=i) } 
#'
#' @import nlme
#' @exportClass gravity
#' @export
gravity <- function(y, x, d, group, data, ln = TRUE, constrained = TRUE, ...) {
  gdata <- data[,c(group, y, x)]
  gdata <- nlme::groupedData( stats::as.formula(paste(paste(y, d, sep=" ~ "), 
                              group, sep=" | ")), data = gdata) 
	  if(ln == TRUE) { 
	   gdata[,x] <- log(abs(gdata[,x]))
	   gdata[,y] <- log(abs(gdata[,y]))
       gdata[gdata == -Inf] <- 0 
       gdata[gdata == Inf] <- 0  		
	  } 		   		
    fmla <- stats::as.formula(paste(paste(y, "~", sep=""), paste(x, collapse= "+")))       	
      if (constrained == FALSE) { 
        print("Running unconstrained gravity model, defaulting to OLS. Please check assumptions")
       gvlmm <- stats::lm(fmla, data = gdata, ...)
	     gvaic <- stats::AIC(gvlmm)
    gm <- list(formula=fmla, gravity=gvlmm, AIC=gvaic, x=gdata[,x], y=gdata[,y],
	           fit = stats::fitted(gvlmm) ) 		   
       } else {	   
	if(!"groupedData" %in% class(gdata)) stop ("Data must be a groupedData object for singly-constrained gravity model")
	  print("Running singly-constrained gravity model")
	    gvlmm <- nlme::lme(fmla, stats::as.formula(paste("random = ~1", group, sep=" | ")), 
		                   data = gdata, ...)
        gvaic <- stats::AIC(gvlmm)
        gm <- list(formula = fmla, gravity = gvlmm, AIC = gvaic, x = gdata[,x], y = gdata[,y], 
	               groups = gdata[,group], fit = stats::fitted(gvlmm) )   	  
    }
	class(gm) <- "gravity"
  return(gm)
}	
