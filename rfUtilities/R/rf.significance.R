#' @title Random Forest model significance test 
#' @description Performs significance test for classification and regression Random Forests models. 
#'
#' @param x randomForest class object
#' @param xdata Independent variables (x) used in model
#' @param p p-value to test for significance in regression models
#' @param q Quantile threshold to test classification models
#' @param nperm Number of permutations
#' @param ... Additional Random Forests arguments 
#'
#' @return A list class object with the following components:
#' @return For Regression problems:    
#' @return     RandRsquare Vector of random R-square values 
#' @return     Rsquare The R-square of the "true" model  
#' @return     Accept Is the model significant at specified p-value (TRUE/FALSE)
#' @return     TestQuantile Quantile threshold used in significance plot
#' @return     pValueThreshold Specified p-value
#' @return     pValue p-values of randomizations
#' @return     nPerm Number of permutations 
#'
#' @return For Classification problems:
#' @return     RandOOB Vector of random out-of-bag (OOB) values 
#' @return     RandMaxError Maximum error of randomizations 
#' @return     test.OOB Error if the "true" model
#' @return     Accept Is the model significant at specified p-value (TRUE/FALSE)
#' @return     TestQuantile Quantile threshold used in significance plot
#' @return     pValueThreshold Specified p-value
#' @return     pValue p-values of randomizations
#' @return     nPerm Number of permutations 
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references
#' Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas connectivity in Yellowstone National Park with landscape genetics. Ecology 91:252-261
#' @references
#' Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer 
#' 
#' @examples 
#' \dontrun{
#' # Regression
#' require(randomForest)
#'   set.seed(1234)	
#'     data(airquality)
#'     airquality <- na.omit(airquality)
#'  ( rf.mdl <- randomForest(x=airquality[,2:6], y=airquality[,1]) )
#'    ( rf.perm <- rf.significance(rf.mdl, airquality[,2:6], nperm=99, ntree=501) )
#'  
#' # Classification
#' require(randomForest)
#'   set.seed(1234)	
#'     data(iris)
#'       iris$Species <- as.factor(iris$Species) 
#'  ( rf.mdl <- randomForest(iris[,1:4], iris[,"Species"], ntree=501) )
#'    ( rf.perm <- rf.significance(rf.mdl, iris[,1:4], nperm=99, ntree=501) ) 
#' }
#'
#' @exportClass significance 
#' @export
rf.significance <- function(x, xdata, q=0.99, p=0.05, nperm=999, ...) 
  {
  if (!inherits(x, "randomForest")) stop("x is not randomForest class object")
    class.p <- function(x, test, nperm) { 
      if ( length( x[x >= test] ) < 1 ) { 
	     error = 1 
	  } else { 
	  	  error = length( x[x < test] ) 
	  }	
      return( error / nperm )
    } 
    regress.p <- function(x, test, nperm) { 
      if ( length( x[x >= test] ) < 1 ) { 
	      error = 1 
	  } else { 
	  	  error = length( x[x >= test] ) 
	 }	
      return( error / nperm )
    } 	
  
  if (x$type == "classification") {  
	test.oob <- stats::median(x$err.rate[,1])
	test.max <- stats::median(max(x$err.rate[,-1]))
	rand.oob <- vector()
	rand.max <- vector()
      for( i in 1:nperm) {	
        rand.y <- sample(x$y, length(x$y)) 
          rf.test <- randomForest::randomForest(x=xdata, y=rand.y, ...)
          rand.oob <- append(rand.oob, stats::median(rf.test$err.rate[,1]) ) 
		rand.max <- append(rand.max, stats::median(max(rf.test$err.rate[,-1])) )
      }
	pValue=round(class.p(x=rand.oob, test=test.oob, nperm=nperm), digits=6)	
    sig <- list(RandOOB=rand.oob, RandMaxError=rand.max, 
                test.OOB=test.oob, test.MaxError=test.max, TestQuantile=q, 
                pValueThreshold=p, pValue=pValue, nPerm=nperm, rf.type=x$type)
      class(sig) <- c("significance","list")			   
	return( sig )				
    }
	
   if (x$type == "regression") {			 
    if (is.factor(x$y)) stop("y cannot be a factor") 
	  test.rsq <- stats::median(x$rsq)
	    test.mse <- stats::median(x$mse)
	  rand.dist <- vector() 
    for( i in 1:nperm) {	
      rand.y <- sample(x$y, length(x$y)) 
        rf.test <- randomForest::randomForest(x=xdata, y=rand.y, ...)
      rand.dist <- append(rand.dist, stats::median(rf.test$rsq)) 
    }	
	pValue=round(regress.p(x=rand.dist, test=test.rsq, nperm=nperm), digits=6)
    sig <- list(RandRsquare=rand.dist, Rsquare=test.rsq, TestQuantile=q, 
                pValueThreshold=p, pValue=pValue, nPerm=nperm, rf.type=x$type) 
	  class(sig) <- c("significance","list")			   
	return( sig )			   
    }
} 
