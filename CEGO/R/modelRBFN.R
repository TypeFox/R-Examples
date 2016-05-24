#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' RBFN Model
#' 
#' Implementation of a distance-based Radial Basis Function Network (RBFN) model, e.g., for mixed or combinatorial input spaces.
#' It is based on employing suitable distance measures for the samples in input space. For reference, see
#' the paper by Moraglio and Kattan (2011).
#'
#' @param x list of samples in input space
#' @param y column vector of observations for each sample
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value, preferably between 0 and 1.
#'      Maximum distances larger 1 are no problem, but may yield scaling bias when different measures are compared.
#' 		Should be non-negative and symmetric. 
#' @param control (list), with the options for the model building procedure:\cr
#' \code{beta} Parameter of the radial basis function: exp(-beta*D), where D is the distance matrix. If beta is not specified, the heuristic in fbeta will be used to determine it, which is default behavior.\cr
#' \code{fbeta} Function f(x) to calculate the beta parameter, x is the maximum distance observed in the input data. Default function is \code{1/(2*(x^2))}. \cr 
#' \code{distances} a distance matrix. If available, this matrix is used for model building, instead of calculating the distance matrix using the parameters \code{distanceFunction}. Default is \code{NULL}.
#'
#'
#' @return a fit (list, modelRBFN), with the options and found parameters for the model which has to be passed to the predictor function:\cr
#' \code{SSQ} Variance of the observations (y)\cr
#' \code{centers} Centers of the RBFN model, samples in input space (see parameters)\cr
#' \code{w} Model parameters (weights) w\cr
#' \code{Phi} Gram matrix \cr
#' \code{Phinv} (Pseudo)-Inverse of Gram matrix\cr
#' \code{w0} Mean of observations (y)\cr
#' \code{dMax} Maximum observed distance\cr
#' \code{D} Matrix of distances between all samples\cr
#' \code{beta} See parameters\cr
#' \code{fbeta} See parameters\cr
#' \code{distanceFunction} See parameters
#' 
#' @seealso \code{\link{predict.modelRBFN}} 
#' 
#' @references Moraglio, Alberto, and Ahmed Kattan. "Geometric generalisation of surrogate model based optimisation to combinatorial spaces." Evolutionary Computation in Combinatorial Optimization. Springer Berlin Heidelberg, 2011. 142-154.
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
#' fit <- modelRBFN(x,y,distancePermutationHamming)
#' #predicted obj. function values
#' ypred <- predict(fit,xtest)$y
#' #plot
#' plot(ytest,ypred,xlab="true value",ylab="predicted value",
#'     pch=20,xlim=c(0.3,1),ylim=c(min(ypred)-0.1,max(ypred)+0.1))
#' abline(0,1,lty=2)
#' @export
###################################################################################
modelRBFN <- function(x,y,distanceFunction,control=list()){ #x sample locations #y observations #distanceFunction function that returns distance matrix
	con<-list( fbeta = function(x) 1/(2*(x^2))
			);
	con[names(control)] <- control;
	control<-con;
	
	#calculate distance matrix?
	if(is.null(control$distances))
		D <-distanceMatrix(x,distanceFunction) 
	else
		D <- control$distances
		
	dMax<-max(D) #maximum distance between samples
	
	if(is.null(control$beta))
		beta <- control$fbeta(dMax) #force a global model, each center has influence everywhere. beta is spread of each RBF. Each now covers the whole search space.
	else
		beta <- control$beta
		
	w0 <- mean(y) #all function values out of reach of any center (see spread due to beta) will be set to the average.

	Phi <- exp(-beta*D^2)

	Phinv <- ginv(Phi) #pseudo inverse (because not guaranteed to be positive definite in non-euclidean space)
	
	w <- Phinv%*%(y-w0) #calculation of wi.		
	
	SSQ <- var(y)
	
	fit <- list(SSQ=SSQ, centers= x, w=w, Phi=Phi, Phinv=Phinv, beta=beta, 
						fbeta=control$fbeta, w0=w0, dMax=dMax, D=D, predAll=FALSE,
						distanceFunction=distanceFunction) #todo not all parameters needed in predictor?
	class(fit) <- "modelRBFN"
	return(fit)
}

###################################################################################
#' Radial Basis Function Network
#'
#' DEPRECATED version of the RBFN model, please use \code{\link{modelRBFN}}
#' 
#' @param x list of samples in input space
#' @param y column vector of observations for each sample
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value
#' @param control options for the model building procedure
#'
#' @keywords internal
#' @export
###################################################################################
combinatorialRBFN <- function(x, y, distanceFunction, control = list()){
	.Deprecated("modelRBFN")
	modelKriging(x,y,distanceFunction,control)
}

###################################################################################
#' Predict: Combinatorial RBFN
#' 
#' Predict with a model fit resulting from \code{\link{modelRBFN}}.
#'
#' @param object fit of the RBFN model (settings and parameters), of class \code{modelRBFN}.
#' @param x list of samples to be predicted
#' @param ... further arguments, not used
#'
#' @return Returned value depends on the setting of \code{object$predAll}\cr
#' TRUE: list with function value (mean) \code{$y} and uncertainty estimate \code{$s} (standard deviation)\cr
#' FALSE:\code{$y}only
#'
#' @seealso \code{\link{modelRBFN}}
#' @export
###################################################################################
predict.modelRBFN <- function(object,x,...){ #x is a new sample, fit is the list of parameters from buildRBFN
	if(!is.list(x))x<-list(x)

	psi <- matrix(unlist(lapply(x,distanceVector,object$centers,object$distanceFunction)),length(object$centers),length(x))
	psi <- exp(-object$beta*psi^2)
	
	pred <-  colSums(apply(psi,2,"*",object$w))+object$w0
	
	res <- list(y=pred)
	if (object$predAll){	
		variance <-  object$SSQ*(1-diag(t(psi)%*%(object$Phinv%*%(psi))))
		s <- sqrt(abs(variance))
		res$s <- as.numeric(s)
	}
	res
}


