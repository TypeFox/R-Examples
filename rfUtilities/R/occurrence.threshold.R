#' @title Test occurrence probability thresholds
#' @description A statistical sensetivity test for occurrence probability thresholds 
#'       
#' @param x             A classification randomForests model object  
#' @param xdata         Independent data used to build model
#' @param class         What class to test
#' @param p             Vector of probability thresholds 
#' @param type          What statistic to use in evaluation ("delta.ss", "sum.ss", "kappa") 
#'
#' @return An "occurrence.threshold" class object contaning a "thresholds" vector object with evaluation statistic and probability thresholds as names.  
#'
#' @details
#' Available threshold evaluation statistics:
#' \itemize{
#' \item   kappa - The Kappa statistic is maximized
#' \item   sum.ss - The sum of sensitivity and specificity is maximized
#' \item   delta.ss - The absolute value of the difference between sensitivity and specificity is minimized
#'  }
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references Jimenez-Valverde, A., & J.M. Lobo (2007). Threshold criteria for conversion of probability of species presence to either-or presence-absence. Acta Oecologica 31(3):361-369 
#' @references Liu, C., P.M. Berry, T.P. Dawson, R.G. Pearson (2005). Selecting thresholds of occurrence in the prediction of species distributions. Ecography 28:385-393.
#'
#' @examples
#' library(randomForest)
#'  data(imports85)
#'   imp85 <- imports85[,-2] 
#'   imp85 <- imp85[complete.cases(imp85), ]
#'   imp85[] <- lapply(imp85, function(x) if (is.factor(x)) x[, drop=TRUE] else x)
#' 
#' y <- ifelse( imp85$numOfDoors != "four", "0", "1")   
#' ( rf.mdl <- randomForest(y = as.factor(y), x = imp85[,-5]) )
#'    ( delta.ss.t <- occurrence.threshold(rf.mdl, imp85[,-5], class = "1") )
#'    ( sum.ss.t <- occurrence.threshold(rf.mdl, imp85[,-5], class = "1", 
#'                                       type = "sum.ss") ) 
#'    ( kappa.ss.t <- occurrence.threshold(rf.mdl, imp85[,-5], class = "1",
#'                                       type = "kappa") )
#'   
#' par(mfrow=c(2,2))
#'   plot(sum.ss.t)
#'   plot(delta.ss.t)
#'   plot(kappa.ss.t)   
#'   
#' @exportClass occurrence.threshold 
#' @export  
occurrence.threshold <- function(x, xdata, class, p = seq(0.10, 0.7, 0.02), 
                                 type = "delta.ss") {
	if(missing(x)) stop( "x (randomForest object) is a required argument") 
  	  if(missing(xdata)) stop( "xdata (data used in rf model) is a required argument") 
	    if(missing(class)) stop( "class is a required argument") 
	if (!inherits(x, "randomForest")) stop("x is not randomForest class object")
	  if (is.null(x$y)) stop("Object does not contain aresponse variable (y) \n")
	    if(x$type != "classification")	stop("Regression not supported \n")	   
          if(length(names(xdata) %in% rownames(x$importance) == TRUE) != length(names(xdata)))
            stop( "varaibles in xdata do not match andom Forests model")	
	          if( (class %in% x$classes) == FALSE) 
                stop("class is  not present in Random Forests model")
	probs <- stats::predict(x, xdata, type="prob")
      probs <- probs[,which(colnames(probs) %in% class)]
	    y <- x$y 
	    y <- ifelse( y == class, 1, 0) 
	test.p <- vector()
      for(i in p) {
          y.p <- ifelse( probs >= i, 1, 0)
        if( type == "kappa") {	
          test.p <- append(test.p, accuracy(y.p, y)$kappa) 
	      }
	    else if (type == "delta.ss") {	
          test.p <- append(test.p, abs(accuracy(y.p, y)$sensitivity - 
	  	                 accuracy(y.p, y)$specificity)) 
	      }
	    else if(type == "sum.ss") {
	      test.p <- append(test.p, (accuracy(y.p, y)$sensitivity + 
	  	                 accuracy(y.p, y)$specificity))
        }
	  }
    names(test.p) <- p
      s <- list(thresholds = test.p, statistic = type)
      class(s) <- c("occurrence.threshold", "list")	  
    return( s )
  }
