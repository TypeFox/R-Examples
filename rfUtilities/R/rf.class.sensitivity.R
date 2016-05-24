#' @title Random Forests class-level sensitivity analysis 
#' @description Performs a sensitivity analysis on a specified class in a random forests model 
#' 
#' @param x  randomForest class object
#' @param xdata  Independent variables used in model
#' @param d  Which class to perturb
#' @param p  Proportion of class to be randomized
#' @param nperm Number of permutations
#' @param plot Plot results (TRUE/FALSE)
#' @param seed  Random seed value
#' @param ...  Additional arguments passed to randomForest  
#'
#' @return List object with following components:
#'   @return mean.error Mean of RMSE
#'   @return sd.error  Standard deviation of RMSE
#'   @return rmse  Root mean squared error (RMSE) for each perturbed probability
#'   @return probs  data.frame with "true" estimate in first column and perturbed probabilities in subsequent columns.        
#'
#' @note Wildlife survey data likely decreases the proportion of imperfect detection (false absences or presences) but can still be a source of error. Because of this it is often necessary to test the model sensitivity of a given class (eg., used verses available habitat). 
#' @note Model sensitivity of false absences is evaluated by randomly assigning a proportion of the specified positive class to the other, refitting the model and estimating the probabilities. Each perturbed estimate is compared against the "true" estimate. Currently only supports binomial models.
#'
#' @author Jeffrey S. Evans  <jeffrey_evans@@tnc.org>
#'
#' @references
#' Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer
#' @references
#' Gardner, R.H., R.V. O'Neill, M.G. Turner, and V.H. Dale (1989). Quantifying scale-dependent effects of animal movements with simple percolation models. Landscape Ecology 3:217-227.
#'
#' @examples
#' library(randomForest)	
#' data(iris)
#'   y <- as.factor(ifelse(iris$Species == "setosa" | 
#'                  iris$Species == "virginica", 1, 0) )
#'     xdata <- iris[,1:4] 
#' 
#' rf.mdl <- randomForest(xdata, y, ntree=501) 
#'   ua <- rf.class.sensitivity(rf.mdl, xdata=xdata, nperm=20, ntree=501, plot=TRUE)
#'       
#' @export        
rf.class.sensitivity <- function(x, xdata, d="1", p=0.05, nperm=999, 
                                 plot=TRUE, seed=NULL, ...) {
  if (!inherits(x, "randomForest")) stop("x is not randomForest class object")
  if (!x$type == "classification") stop( "x is not a classification object")
  if(!is.null(seed)) { set.seed(seed) } else { set.seed(.Random.seed[1]) }
    rmse <- function(o,p) sqrt( mean( (o - p )^2 ) )	
	values <- as.character(unique(x$y))
	pred <- stats::predict(x, xdata, type="prob")
	nc <- which(colnames(pred) == d)
	mpred <- data.frame(obs=pred[,nc])
	names(mpred) <- "obs" 
	if(plot == TRUE) 
	  graphics::plot(stats::density(mpred[,1]), main="Perturbed model probabilities", type="n")
	    for(i in 1:nperm) {
	      y <- x$y
            samp <- sample(which( y %in% d ), length(which( y %in% d )) * p)
	        y[samp] <- values[values != d] 
	    	pmdl <- randomForest::randomForest(xdata, y)     
	    	mpred <- data.frame(mpred, stats::predict(pmdl, xdata, type="prob")[,nc])			  
	    	  if(plot == TRUE) {
			    pden <- stats::density(mpred[,i+1])
			    graphics::lines(pden$x, pden$y, col="grey")
			  }
	    }
	      if(plot == TRUE) {
            d <- stats::density(mpred[,1])
	        graphics::lines(d$x, d$y, col="red", lwd=2)
	        graphics::legend("topleft", legend=c("observed", "perturbed"), 
	                         col=c("red", "grey"), lwd=c(2,1), bg="white")
          }		
            names(mpred)[2:ncol(mpred)] <- paste0("sim", seq(1:(ncol(mpred)-1))) 
    	    error <- vector()
        for(i in 2:ncol(mpred)) error <- append(error, rmse(mpred[,1],mpred[,i]))
      cat("Mean error: ", mean(error), "\n")	
    cat("Standard deviation of error: ", stats::sd(error), "\n")	 
  list(mean.error = mean(error), sd.error = stats::sd(error), rmse = error, probs = mpred)
  }
