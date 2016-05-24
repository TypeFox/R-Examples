#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Distance based Linear Model
#' 
#' A simple linear model based on arbitrary distances. Comparable to a k nearest neighbor model, but potentially able to extrapolate
#' into regions of improvement. Used as a simple baseline by Zaefferer et al.(2014).
#'
#' @param x list of samples in input space
#' @param y  matrix, vector of observations for each sample
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value, preferably between 0 and 1.
#'      Maximum distances larger 1 are no problem, but may yield scaling bias when different measures are compared.
#' 		Should be non-negative and symmetric. 
#' @param control currently unused, defaults to \code{list()}
#'
#'
#' @return a fit (list, modelLinear), with the options and found parameters for the model which has to be passed to the predictor function:\cr
#' \code{x} samples in input space (see parameters)\cr
#' \code{y} observations for each sample (see parameters)\cr
#' \code{distanceFunction} distance function (see parameters)
#' 
#' @seealso \code{\link{predict.modelLinear}} 
#' 
#' @references Zaefferer, Martin; Stork, Joerg; Friese, Martina; Fischbach, Andreas; Naujoks, Boris; Bartz-Beielstein, Thomas. (2014). Efficient global optimization for combinatorial problems. In Proceedings of the 2014 conference on Genetic and evolutionary computation (GECCO '14). ACM, New York, NY, USA, 871-878. DOI=10.1145/2576768.2598282 http://doi.acm.org/10.1145/2576768.2598282 
#'
#' @examples
#' #set random number generator seed
#' set.seed(1)
#' #simple test landscape
#' fn <- landscapeGeneratorUNI(1:5,distancePermutationHamming)
#' #generate data for training and test
#' x <- unique(replicate(40,sample(5),FALSE))
#' xtest <- x[-(1:15)]
#' x <- x[1:15]
#' #determin true objective function values
#' y <- fn(x)
#' ytest <- fn(xtest)
#' #build model
#' fit <- modelLinear(x,y,distancePermutationHamming)
#' #predicted obj. function values
#' ypred <- predict(fit,xtest)$y
#' #plot
#' plot(ytest,ypred,xlab="true value",ylab="predicted value",
#'     pch=20,xlim=c(0.3,1),ylim=c(min(ypred)-0.1,max(ypred)+0.1))
#' abline(0,1,lty=2)
#' @export
###################################################################################
modelLinear <- function(x, y, distanceFunction, control=list()){ #linear distance model
	fit<-list()
	fit$distanceFunction <- distanceFunction
	fit$x <- x
	fit$y <- y
	class(fit) <- "modelLinear"
	fit
}

###################################################################################
#' Linear Distance-Based Model
#'
#' DEPRECATED version of the linear, distance-based model, please use \code{\link{modelLinear}}
#' 
#' @param x list of samples in input space
#' @param y column vector of observations for each sample
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value
#' @param control options for the model building procedure
#'
#' @keywords internal
#' @export
###################################################################################
combinatorialLM <- function(x, y, distanceFunction, control = list()){
	.Deprecated("modelLinear")
	modelKriging(x,y,distanceFunction,control)
}

###################################################################################
#' Predict: Combinatorial Kriging
#' 
#' Predict with amodelLinear fit.
#'
#' @param object fit of the Kriging model (settings and parameters), of class \code{modelLinear}.
#' @param x list of samples to be predicted
#' @param ... further arguments, not used
#'
#' @return numeric vector of predictions
#'
#' @seealso \code{\link{modelLinear}}
#' @export
###################################################################################
predict.modelLinear <- function(object,x,...){ #approach: sort by distance, build linear model, predict at zero.
	if(!is.list(x))x<-list(x)
	pred=NULL
	for(i in 1:length(x)){
		distx <- distanceVector(x[[i]],object$x,object$distanceFunction)

		fy=object$y
		dx=distx #global model

		s2 <- sort(unique(distx))[2]
		index <- which(distx <= s2)
		fy=fy[index]
		dx=dx[index] #local model

		index2 <- which(dx < s2)
		fy1 = mean(fy[index2])
		fy2 = mean(fy[-index2])
		difY = fy2-fy1
		difD = max(dx)-min(dx)
		m=difY/difD
		pred <- c(pred, fy1-min(dx)*m)

	}
	list(y=pred,s=rep(0,length(pred)))
}