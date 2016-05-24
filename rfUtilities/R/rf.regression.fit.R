#' @title Random Forest fit statistics  
#' @description Evaluatues fit and overfit of random forests regression 
#' 
#' @param x    randomForest regression object
#' @return     A list and rf.fit class object with "fit" matrix of fit statistics and "message" indicating overfit risk. 
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @examples
#'   library(randomForest)
#'   set.seed(131)
#'   data(airquality)
#'   airquality <- na.omit(airquality)
#'   ( rf.aq <- randomForest(airquality[,1:3], airquality[,"Ozone"]) )	
#'   rf.regression.fit(rf.aq)
#'
#' @export
rf.regression.fit <- function(x) {
  if(!x$type == "regression") stop("Classification models not supported")
  rmse <- sqrt(mean((x$predicted - x$y)^2))        
  r2 <- mean(x$rsq)                                
  f2 <- r2 / (1 - r2)                              
  n <- length(x$y)
  k <- nrow(x$importance)  
  overfit.ratio <- round(n / k)  # overfitting ratio, ideal is > 10
  fit.Actuals.pred <- cbind(x$predicted, x$y)
  accuracy <- stats::median(apply(fit.Actuals.pred, 1, min) / 
                            apply(fit.Actuals.pred, 1, max))
    if(overfit.ratio > 9) { 
	  msg <- "Model is not overfit" 
       } else { 
	  msg <- "Model may be overfit" 
	}
	fit.stats <- t(data.frame(RMSE = round(rmse,3), 
	                          R.squared = round(r2,3), 
						      Cohen.f2 = round(f2,3),
                     	      Accuracy = round(accuracy,3), 
						      Overfitting.ratio = round(overfit.ratio,3)))
	fit <- list(fit = fit.stats, message = msg)
    class(fit) <- c("rf.fit","list")	
  return(fit)                         
}
