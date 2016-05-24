#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Unimodal Fitness Landscape
#' 
#' This function generates uni-modal fitness landscapes based on distance measures.
#' The fitness is the distance to a reference individual. Hence, the reference individual
#' is the optimum of the landscape. 
#'
#' @param ref Reference individual
#' @param distanceFunction Distance function, used to evaluate d(x,ref), where x is an arbitrary new individual
#'
#' @return returns a function. The function requires a list of candidate solutions as its input, where each solution is suitable for use with the distance function. The function returns a numeric vector.
#'
#' @references Moraglio, Alberto, Yong-Hyuk Kim, and Yourim Yoon. "Geometric surrogate-based optimisation for permutation-based problems." Proceedings of the 13th annual conference companion on Genetic and evolutionary computation. ACM, 2011.
#'
#' @seealso \code{\link{landscapeGeneratorMUL}}, \code{\link{landscapeGeneratorGaussian}} 
#'
#' @examples
#' fun <- landscapeGeneratorUNI(ref=1:7,distancePermutationCos)
#' ## for single solutions, note that the function still requires list input:
#' x <- 1:7
#' fun(list(x))
#' x <- 7:1
#' fun(list(x))
#' x <- sample(7)
#' fun(list(x))
#' ## multiple solutions at once:
#' x <- replicate(5,sample(7),FALSE)
#' fun(x)
#'
#' @export
###################################################################################
landscapeGeneratorUNI <- function(ref=1:10,distanceFunction){  #N= number of elements of a string, K= number of neighbours, PI = relative position of each neighbour in string, g= set of fitness functions for each combination of string components
	distanceFunction #lazy evaluation fix, faster than force()
	ref #lazy evaluation fix, faster than force()
	function(x){
		distanceVector(ref,x,distanceFunction)
	}
}


###################################################################################
#' Multimodal Fitness Landscape
#' 
#' This function generates multi-modal fitness landscapes based on distance measures.
#' The fitness is the minimal distance to several reference individuals. Hence, each reference individual
#' is an optimum of the landscape. 
#'
#' @param ref list of reference individuals
#' @param distanceFunction Distance function, used to evaluate d(x,ref[[n]]), where x is an arbitrary new individual
#'
#' @return returns a function. The function requires a list of candidate solutions as its input, where each solution is suitable for use with the distance function. The function returns a numeric vector.
#'
#' @seealso \code{\link{landscapeGeneratorUNI}}, \code{\link{landscapeGeneratorGaussian}}
#'
#' @examples
#' fun <- landscapeGeneratorMUL(ref=list(1:7,c(2,4,1,5,3,7,6)),distancePermutationCos)
#' x <- 1:7
#' fun(list(x))
#' x <- c(2,4,1,5,3,7,6)
#' fun(list(x))
#' x <- 7:1
#' fun(list(x))
#' x <- sample(7)
#' fun(list(x))
#' ## multiple solutions at once:
#' x <- append(list(1:7,c(2,4,1,5,3,7,6)),replicate(5,sample(7),FALSE))
#' fun(x)
#'
#' @export
###################################################################################
landscapeGeneratorMUL <- function(ref=list(1:10),distanceFunction){
	distanceFunction #lazy evaluation fix, faster than force()
	ref #lazy evaluation fix, faster than force()
	function(x){
		k = length(ref)
		d <- matrix(NA,k,length(x))#rep(NA,k)
		for(i in 1:k)
			d[i,] <- distanceVector(ref[[i]],x,distanceFunction)
		do.call(pmin.int, lapply(1:nrow(d), function(i)d[i,])) #fast row minimum for matrix d, the return value for each candidate solution
	}	
}


###################################################################################
#' Create Gaussian Landscape
#' 
#' This function is loosely based on the Gaussian Landscape Generator by Bo Yuan and Marcus Gallagher.
#' It creates a Gaussian Landscape every time it is called. This Landscape can be evaluated like a function.
#' To adapt to combinatorial spaces, the Gaussians are here based on a user-specified distance measure.
#' Due to the expected nature of combinatorial spaces and their lack of direction, the resulting
#' Gaussians are much simplified in comparison to the continuous, vector-valued case (e.g., no rotation).
#' Since the \code{CEGO} package is tailored to minimization, the landscape is inverted.
#'
#' @param nGaussian number of Gaussian components in the landscape. Default is 10.
#' @param upper upper boundary of the \code{distanceFunction}. Controls width of Gaussian components. Default is 1.
#' @param ratio minimal function value of the local minima. Default is 0.2. (Note: Global minimum will be at zero, local minimal will be in range \code{[ratio;1]})
#' @param seed seed for the random number generator used before creation of the landscape. Generator status will be saved and reset afterwards.
#' @param distanceFunction A function of type \code{f(x,y)}, to evaluate distance between to samples in their given representation.
#' @param creationFunction function to randomly generate the centers of the Gaussians, in form of their given representation.
#'
#' @return returns a function.The function requires a list of candidate solutions as its input, where each solution is suitable for use with the distance function.
#'
#' @references B. Yuan and M. Gallagher (2003) "On Building a Principled Framework for Evaluating and Testing Evolutionary Algorithms: A Continuous Landscape Generator". 
#' In Proceedings of the 2003 Congress on Evolutionary Computation, IEEE, pp. 451-458, Canberra, Australia.
#'
#' @examples 
#' #rng seed
#' seed=101
#' # distance function
#' dF <- function(x,y)(sum((x-y)^2)) #sum of squares 
#' #dF <- function(x,y)sqrt(sum((x-y)^2)) #euclidean distance
#' # creation function
#' cF <- function()runif(1)
#' # plot pars
#' par(mfrow=c(3,1),mar=c(3.5,3.5,0.2,0.2),mgp=c(2,1,0))
#' ## uni modal distance landscape
#' # set seed
#' set.seed(seed)
#' #landscape
#' lF <- landscapeGeneratorUNI(cF(),dF)
#' x <- seq(from=-0,by=0.001,to=1)
#' plot(x,lF(x),type="l")
#' ## multi-modal distance landscape
#' # set seed
#' set.seed(seed)
#' #landscape
#' lF <- landscapeGeneratorMUL(replicate(5,cF(),FALSE),dF)
#' x <- seq(from=-0,by=0.001,to=1)
#' plot(x,lF(x),type="l")
#' ## glg landscape
#' #landscape
#' lF <- landscapeGeneratorGaussian(nGaussian=20,upper=1,
#' ratio=0.8,seed=seed,dF,cF)
#' x <- seq(from=-0,by=0.001,to=1)
#' plot(x,lF(x),type="l")
#'
#' @export
###################################################################################
landscapeGeneratorGaussian <- function(nGaussian=10,upper=1,ratio=0.2,seed=1, distanceFunction, creationFunction){
	## save seed status
	if(exists(as.character(substitute(.Random.seed))))
		SAVESEED<-.Random.seed
	else
		SAVESEED=NULL
		
	## set seed
	set.seed(seed)
	
	## create landscape
	fit <- landscapeGeneratorGaussianBuild(nGaussian,upper,ratio,creationFunction)
	fit$df <- distanceFunction 
	
	## load seed status
	if(!is.null(SAVESEED))
		assign(".Random.seed", SAVESEED, envir=globalenv())
		
	fit	
	## create output function	
	fun <- function(x){
		landscapeGeneratorGaussianEval(x,fit)$value
	}
	attributes(fun) <- fit
	fun
}



###################################################################################
#' Gaussian Landscape Core function
#' 
#' Core Gaussian landscape function. Should not be called directly, as it does not contain proper seed handling.
#'
#' @param nGaussian number of Gaussian components in the landscape. Default is 10.
#' @param upper upper limit of the distance function to be used, controls the width of Gaussians (variance). Default is 1.
#' @param ratio minimal function value of the local minima. Default is 0.2. (Note: Global minimum will be at zero, local minimal will be in range \code{[ratio;1]})
#' @param creationFunction function to randomly generate the centers of the gaussians
#'
#' @return returns a list, with the following items:\cr
#' \code{centers} samples which are the centers of each Gaussian
#' \code{covinv} inverse of variance of each Gaussian 
#' \code{opt} value at randomly chosen optimum center
#' \code{nGauss} number of Gaussian components
#'
#' @keywords internal
#' @export
#' @seealso \code{\link{landscapeGeneratorGaussian}} 
###################################################################################
landscapeGeneratorGaussianBuild <- function(nGaussian=10,upper=1, ratio=0.2,creationFunction){
	ratio <- 1-ratio
	if (nGaussian<=0|ratio<=0|ratio>=1){
		stop('Incorrect parameter values for gaussian landscape generator')
	}

	varianceRange=upper/20 
	
	#variance=matrix(runif(nGaussian*dimension),nGaussian,dimension)*variancerange+0.05;  # add 0.05 to avoid zero variance     #TODO or mat mult?
	variance <- matrix(runif(nGaussian),nGaussian,1)*varianceRange+0.05;  # add 0.05 to avoid zero variance    

	covmatrix <- variance 
	covmatrix_inv <- 1/covmatrix
	
	# Generate randomly the centers of the gaussians
	meanVector <- replicate(nGaussian,creationFunction(),simplify=FALSE)

	# assign values to components
	optimumValue=rep(0,nGaussian) #initialize
	optimumValue[1]=1     # the first Gaussian is set to be the global optimum

	# values of others are randomly generated within [0,ratio]
	optimumValue[2:nGaussian]=matrix(runif(1*(nGaussian-1)),1,nGaussian-1)*ratio
	list(centers=meanVector,covinv=covmatrix_inv,opt=1-optimumValue,nGauss=nGaussian)
}

###################################################################################
#' Gaussian Landscape Evaluation
#' 
#' Evaluate a Gaussian landscape. Should not be called directly.
#'
#' @param x list of samples to evaluate 
#' @param glg list of values defining the Gaussian Landscape, created by \code{landscapeGeneratorGaussianBuild}.
#'
#' @return returns a list, with the following items:\cr
#' \code{value} value of the combined landscape
#' \code{components} value of each component
#'
#' @keywords internal
#' @export
#' @seealso \code{\link{landscapeGeneratorGaussian}} 
###################################################################################
landscapeGeneratorGaussianEval <- function(x,glg){
	covmatrix_inv <- glg$covinv #the inverse covariance matrix of each component
	meanvector <- glg$centers     #centers of each Gaussian component
	optimumvalue <- 1-glg$opt   #the peak value of each component
	nGaussian <-  glg$nGauss  #total number of components
	
	p <- length(x)                 # p: number of individuals;   
	#if(is.null(p))p<-1	
	tmp <- matrix(0,nGaussian,p)

	#----------------------------------------------------
	for(i in 1:nGaussian) {             # calculate the values generated by each component		
		newx <- distanceVector(meanvector[[i]],x,glg$df)		#x-t(matrix(meanvector[i,],length(meanvector[i,]),p,byrow=FALSE))
		tmp[i,] <- covmatrix_inv[i]*newx     	   
	}
	f <- exp(-0.5*tmp)            # f is a nGaussian-by-p matrix

	f <- f*matrix(optimumvalue,length(optimumvalue),p,byrow=FALSE)# multiply the peak value of each component
	                # the value of each individual generated by each component
	value <- apply(f,2,max)   # choose the maximum values as the fitness values
	list(value=1-value,components=f)
}
