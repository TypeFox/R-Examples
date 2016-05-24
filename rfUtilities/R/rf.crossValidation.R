#' @title Random Forest Classification Model Cross-validation 
#' @description Implements a permutation test cross-validation for Random Forests classification models
#'    
#' @param x                 random forest object
#' @param xdata             x data used in model
#' @param p                 Percent data withhold
#' @param n                 Number of cross validations
#' @param seed              Sets random seed in R global environment
#' @param ...               Additional arguments passed to Random Forests 
#'
#' @return  A "rf.cv" class object with the following components:
#' \itemize{ 
#' \item  cross.validation$cv.users.accuracy        Class-level users accuracy for the subset cross validation data   
#' \item  cross.validation$cv.producers.accuracy    Class-level producers accuracy for the subset cross validation data    
#' \item  cross.validation$cv.oob                   Global and class-level OOB error for the subset cross validation data    
#' \item  model$model.users.accuracy                Class-level users accuracy for the model 
#' \item  model$model.producers.accuracy            Class-level producers accuracy for the model 
#' \item  model$model.oob                           Global and class-level OOB error for the model  
#'  }
#'
#' @details
#' The crossvalidation statistics are based on the prediction error on the witheld data: 
#' 
#' Total observed accuracy represents the percent correctly classified (AKA, ) and is considered as a naive measure of agreement. The diagonal of the confusion matrix represents correctly classified observations where off-diagonals represent cross-classification error. The primary issue with this evaluation is that does not reveal if error was evenly distributed between classes.   
#'
#' To represent the balance of error one can use omission and commission statistics such as estimates of users and producers accuracy. User's accuracy corresponds to error of commission (inclusion), observations being erroneously included in a given class. The commission errors are represented by row sums of the matrix. Producer's accuracy corresponds to error of omission (exclusion), observations being erroneously excluded from a given class. The omission errors are represented by column sums of the matrix.
#'
#' None of the previous statistics account for random agreement influencing the accuracy measure. The kappa statistic is a chance corrected metric that reflects the difference between observed agreement and agreement expected by random chance. A kappa of k=0.85 would indicate that there is 85% better agreement than by chance alone.  
#'
#' \itemize{ 
#' \item   pcc = [Number of correct observations / total number of observations] 
#' \item   pcc = [Number of correct observations / total number of observations] 
#' \item   producers accuracy =  [Number of correct / total number of correct and omission errors] 
#' \item   k = (observed accuracy - chance agreement) / (1 - chance agreement) where; change agreement = sum[product of row and column totals for each class] 
#'  }
#'
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references
#' Evans, J.S. and S.A. Cushman (2009) Gradient Modeling of Conifer Species Using Random Forest. Landscape Ecology 5:673-683.
#' @references
#' Murphy M.A., J.S. Evans, and A.S. Storfer (2010) Quantify Bufo boreas connectivity in Yellowstone National Park with landscape genetics. Ecology 91:252-261
#' @references
#' Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer 
#' 
#' @examples 
#' require(randomForest)
#'   data(iris)
#'     iris$Species <- as.factor(iris$Species)    	
#'       set.seed(1234)	
#' ( rf.mdl <- randomForest(iris[,1:4], iris[,"Species"], ntree=501) )
#'   ( rf.cv <- rf.crossValidation(rf.mdl, iris[,1:4], p=0.10, n=99, ntree=501) )
#'
#'    # Plot cross validation verses model producers accuracy
#'    par(mfrow=c(1,2)) 
#'      plot(rf.cv, type = "cv", main = "CV producers accuracy")
#'      plot(rf.cv, type = "model", main = "Model producers accuracy")
#'
#'    # Plot cross validation verses model oob
#'    par(mfrow=c(1,2)) 
#'      plot(rf.cv, type = "cv", stat = "oob", main = "CV oob error")
#'      plot(rf.cv, type = "model", stat = "oob", main = "Model oob error")	  
#'	  
#' @exportClass rf.cv
#' @export 	
rf.crossValidation <- function(x, xdata, p=0.10, n=99, seed=NULL, ...) {
  if (!inherits(x, "randomForest")) stop("x is not randomForest class object")
    if (!x$type == "classification") stop( "Random Forests model is not a classification") 
      if(!is.null(seed)) { set.seed(seed) }	
    classes <- as.vector(levels( x$y ))	
      sample.size = round( (length(x$y) * p) / length(x$classes), digits=0) 
	    cv.ua <- as.data.frame(array(0, dim=c(0,length(classes))))
	      cv.pa <- as.data.frame(array(0, dim=c(0,length(classes))))
	      mdl.ua <- as.data.frame(array(0, dim=c(0,length(classes)))) 
	    mdl.pa <- as.data.frame(array(0, dim=c(0,length(classes))))
      mdl.oob <- as.data.frame(array(0, dim=c(0,2+length(classes))))
	cv.oob <- as.data.frame(array(0, dim=c(0,2+length(classes))))		
      for(i in 1:n) {	
	    samp.index <- vector()
          for(r in unique(x$y)) {			  
            cvalue <- which(x$y == r) 	  
            samp.index <- append(samp.index, sample(cvalue, sample.size)) 
          }		  
	      tx <- xdata[samp.index,]
	      ty <- x$y[samp.index]
	      mx <- xdata[-samp.index,]
	      my <- x$y[-samp.index] 
          rf.fit <- randomForest::randomForest(y=as.factor(my), x=mx, ytest=as.factor(ty), xtest=tx, ...)    
        
        cv.acc <- accuracy(rf.fit$test$predicted,  ty)
          if(!length(classes) == length(unique(ty))) {
	        #### add code to check for presence of classes and set to NA if absent classes occur   
		  }		  
		cv.ua <- rbind(cv.ua, cv.acc$users.accuracy) 
		  cv.pa <- rbind(cv.pa, cv.acc$producers.accuracy) 
	        cv.oob <- rbind(mdl.oob, c(apply(rf.fit$test$err.rate, MARGIN = 2, stats::median),
		                    cv.acc$kappa))	  
		
		mdl.acc <- accuracy(rf.fit$predicted, my) 
		  if(!length(classes) == length(unique(my))) {
	        #### add code to check for presence of classes and set to NA if absent classes occur   
		  }
	    mdl.ua <- rbind(mdl.ua, mdl.acc$users.accuracy) 
		  mdl.pa <- rbind(mdl.pa, mdl.acc$producers.accuracy) 
            mdl.oob <- rbind(mdl.oob, c(apply(rf.fit$err.rate, MARGIN = 2, stats::median),
		                     mdl.acc$kappa))
	}  
    names(cv.ua) <- c(classes) 
      names(cv.pa)  <- c(classes) 
        names(mdl.ua) <- c(classes) 	  
          names(mdl.pa) <- c(classes) 	  
            names(mdl.oob) <- c("OOB", classes, "kappa") 
              names(cv.oob) <- c("OOB", classes, "kappa")	  
    acc <- list( cross.validation = list(cv.users.accuracy = cv.ua, 
	            cv.producers.accuracy=cv.pa, cv.oob = cv.oob), 
	            model = list(model.users.accuracy = mdl.ua, 
				             model.producers.accuracy=mdl.pa, 
			    			 model.oob=mdl.oob) )  	
    class(acc) <- c("rf.cv", "list")
  return( acc )
} 
