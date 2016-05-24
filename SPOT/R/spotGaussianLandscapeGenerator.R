###################################################################################
#' Create Gaussian Landscape (multiple)
#' 
#' This function is based on the Gaussian Landscape Generator by Bo Yuan and Marcus Gallagher.
#' The difference to \code{\link{spotGlgCreate}} is, that this function creates multiple landscapes at once
#' which can be used as a multi objective test problem. They are simply created sequentially after the random
#' number generator seed is set. Any differences are due to their non-deterministic nature.
#'
#' @param oDimension dimensionality of output, i.e. number of landscapes. Default is 2.
#' @param iDimension dimensionality of the landscapes input space. Default is 2.
#' @param nGaussian number of Gaussian components in the landscape. Default is 10.
#' @param lower lower boundary of the landscape, defaults to \code{rep(0,dimension)}.
#' @param upper upper boundary of the landscape, defaults to \code{rep(1,dimension)}.
#' @param globalvalue the global maximum value, i.e. the maximum of the Gaussian component with the largest value. Default is 1.
#' @param ratio maximum ratio of the local maxima, local optima are randomly generated within \code{[0,globalvalue*ratio]}. Has to be larger than 0 and smaller than 1. Defaults to 0.8.
#' @param seed seed for the random number generator used before creation of the landscape. Generator status will be saved and reset afterwards.
#'
#' @return returns a function. The function takes a point (vector) as input, with as many values as specified in \code{dimension}.
#' The function returns a vector value, representing the output from each landscape.
#'
#' @references B. Yuan and M. Gallagher (2003) "On Building a Principled Framework for Evaluating and Testing Evolutionary Algorithms: A Continuous Landscape Generator". 
#' In Proceedings of the 2003 Congress on Evolutionary Computation, IEEE, pp. 451-458, Canberra, Australia.
#'
#' @seealso \code{\link{spotGlgCreate}}, \code{\link{spotGlgCreateRot}}, \code{\link{spotGlgCreateRotSearched}} 
#'
#' @export
###################################################################################
spotGlgCreateN<- function(oDimension=2,iDimension=2,nGaussian=rep(10,oDimension),
				lower=rep(0,iDimension),upper=rep(1,iDimension),
				globalvalue=rep(1,oDimension),ratio=rep(0.8,oDimension),seed=1){
	fit <- spotGlgInitN(oDimension,iDimension,nGaussian,lower,upper,globalvalue,ratio,seed)
	fun <- function(x){
		spotGlgEvalN(x,fit)$value
	}
	#attributes(fun) <- fit
	fun
}

###################################################################################
#' Initialize Gaussian Landscape (multiple)
#' 
#' This function is based on the Gaussian Landscape Generator by Bo Yuan and Marcus Gallagher.
#' It randomly initializes multiple Gaussian Landscapes with the specified parameters.
#' Does not have to be called by user, if \code{\link{spotGlgCreateN}} is used.
#'
#' @param oDimension dimensionality of output, i.e. number of landscapes. Default is 2.
#' @param iDimension dimensionality of the landscapes input space. Default is 2.
#' @param nGaussian number of Gaussian components in the landscape. Default is 10.
#' @param lower lower boundary of the landscape, defaults to \code{rep(0,dimension)}.
#' @param upper upper boundary of the landscape, defaults to \code{rep(1,dimension)}.
#' @param globalvalue the global maximum value, i.e. the maximum of the Gaussian component with the largest value. Default is 1.
#' @param ratio maximum ratio of the local maxima, local optima are randomly generated within \code{[0,globalvalue*ratio]}. Has to be larger than 0 and smaller than 1. Defaults to 0.8.
#' @param seed seed for the random number generator used before creation of the landscape. Generator status will be saved and reset afterwards.
#'
#' @return returns a list containing again lists for each landscape,
#'  containing with the following items:\cr
#' \code{mean} Matrix containing the mean vectors of the Gaussian components in the landscape, i.e. the locations of the local maxima of the functions. First vector (i.e. first row) will be the global maximum.
#' \code{covinv} Inverse of covariance matrix of each Gaussian component, stored as 3-dimensional array.
#' \code{opt} optimal values, i.e. maxima of the Gaussian components
#' \code{ngauss} number of Gaussian components
#' \code{d} is the \code{dimension}
#'
#' @references B. Yuan and M. Gallagher (2003) "On Building a Principled Framework for Evaluating and Testing Evolutionary Algorithms: A Continuous Landscape Generator". 
#' In Proceedings of the 2003 Congress on Evolutionary Computation, IEEE, pp. 451-458, Canberra, Australia.
#'
#' @seealso \code{\link{spotGlgCreateN}}, \code{\link{spotGlgEvalN}}, \code{\link{spotGlgInit}}
#'
#' @export
###################################################################################
spotGlgInitN <- function(oDimension=2,iDimension=2,nGaussian=rep(10,oDimension),
				lower=rep(0,iDimension),upper=rep(1,iDimension),
				globalvalue=rep(1,oDimension),ratio=rep(0.8,oDimension),seed=1){
	if(exists(as.character(substitute(.Random.seed))))
		SAVESEED<-.Random.seed
	else
		SAVESEED=NULL
	set.seed(seed)
	
	res <- list()
	for(i in 1:oDimension){
		res[[i]] <- spotGlgCore(iDimension,nGaussian[i],lower,upper,globalvalue[i],ratio[i])
	}
	
	if(!is.null(SAVESEED))
		assign(".Random.seed", SAVESEED, envir=globalenv())
	res
}

###################################################################################
#' Gaussian Landscape Evaluation (multiple)
#' 
#' This function is based on the Gaussian Landscape Generator by Bo Yuan and Marcus Gallagher.
#' It randomly evaluates one or several points in several Gaussian Landscapes at once, as created by \code{spotGlgInitN}.
#'
#' @param x matrix of sample sites, containing one point in each row.
#' @param glglist list of lists (one for each landscape) of values defining the Gaussian Landscape, created by \code{spotGlgInit}.
#'
#' @return returns a list, with the following items:\cr
#' \code{value} vector of values from each landscape.
#' \code{components} value of each component (will be NA)
#'
#' @references B. Yuan and M. Gallagher (2003) "On Building a Principled Framework for Evaluating and Testing Evolutionary Algorithms: A Continuous Landscape Generator". 
#' In Proceedings of the 2003 Congress on Evolutionary Computation, IEEE, pp. 451-458, Canberra, Australia.
#'
#' @seealso \code{\link{spotGlgCreateN}}, \code{\link{spotGlgInitN}}, \code{\link{spotGlgEval}} , \code{\link{spotGlgEvalRot}}  
#'
#' @export
###################################################################################
spotGlgEvalN <- function(x,glglist){
	value <- NULL
	for(i in 1:length(glglist)){
		res <- spotGlgEval(x,glglist[[i]])
		value <- cbind(value, res$value)
	}
	list(value=value,components=NA)
}

###################################################################################
#' Create Gaussian Landscape (rotated) with random search
#' 
#' This function is very similar to \code{\link{spotGlgCreateRot}}.
#' However, to allow for meaningfull results when comparing
#' hypervolumes, this function supports random-searching the created landscapes.
#' The resulting mean hypervolume of the random-search can be used for comparisons.
#'
#' @param iDimension dimensionality of the landscapes input space. Default is 2.
#' @param alpha rotation angle. Default is pi/2.
#' @param nGaussian number of Gaussian components in the landscape. Default is 10.
#' @param lower lower boundary of the landscape, defaults to \code{rep(0,dimension)}.
#' @param upper upper boundary of the landscape, defaults to \code{rep(1,dimension)}.
#' @param globalvalue the global maximum value, i.e. the maximum of the Gaussian component with the largest value. Default is 1.
#' @param ratio maximum ratio of the local maxima, local optima are randomly generated within \code{[0,globalvalue*ratio]}. Has to be larger than 0 and smaller than 1. Defaults to 0.8.
#' @param seed seed for the random number generator used before creation of the landscape. Generator status will be saved and reset afterwards.
#' @param repeats number of random-search runs performed on the landscape.
#' @param evals evaluations of the landscapes in each run of the random-search.
#'
#' @return returns a function. The function takes a point (vector) as input, with as many values as specified in \code{dimension}.
#' The function returns a vector value, representing the output from the two landscapes landscape. The mean hypervolume
#' is the "meanhvol" attribute of said function.
#'
#' @seealso \code{\link{spotGlgCreate}} \code{\link{spotGlgCreateN}} \code{\link{spotGlgCreateRot}}
#'
#' @export
###################################################################################
spotGlgCreateRotSearched <-function(iDimension=2,alpha=pi/2,nGaussian=10,
				lower=rep(0,iDimension),upper=rep(1,iDimension),
				globalvalue=1,ratio=0.8,seed=1,repeats=100,evals=1000){
	####
	##set seed and save rng status
	if(exists(as.character(substitute(.Random.seed))))
		SAVESEED<-.Random.seed
	else
		SAVESEED=NULL
	set.seed(seed)	
	####
	## create function
	fn <- spotGlgCreateRot(iDimension,alpha,nGaussian,lower,upper,globalvalue,ratio,seed)
	fun <- function(x) globalvalue -fn(x)
	attributes(fun)=attributes(fn)
	####
	##calculate random search hyper volume
	hvol = NULL
	for(i in 1:repeats){
		xran = matrix(runif(iDimension*evals),,evals)
		xran = t(xran*(upper-lower)+lower)
		yran = fun(xran)
		hvol = c(hvol,dominated_hypervolume(t(yran),ref=c(globalvalue,globalvalue)))
	}
	avghvol=mean(hvol)
	attr(fun,"meanhvol")=avghvol
	####
	##reset rng
	if(!is.null(SAVESEED))
		assign(".Random.seed", SAVESEED, envir=globalenv())	
	fun
}

###################################################################################
#' Create Gaussian Landscape (rotated)
#' 
#' This function is based on the Gaussian Landscape Generator by Bo Yuan and Marcus Gallagher.
#' The difference to \code{\link{spotGlgCreate}} is, that this function creates one landscape
#' and a second one is created by rotating the first. This can be used as a simple two-objective
#' test function.
#'
#' @param iDimension dimensionality of the landscapes input space. Default is 2.
#' @param alpha rotation angle. Default is pi/2.
#' @param nGaussian number of Gaussian components in the landscape. Default is 10.
#' @param lower lower boundary of the landscape, defaults to \code{rep(0,dimension)}.
#' @param upper upper boundary of the landscape, defaults to \code{rep(1,dimension)}.
#' @param globalvalue the global maximum value, i.e. the maximum of the Gaussian component with the largest value. Default is 1.
#' @param ratio maximum ratio of the local maxima, local optima are randomly generated within \code{[0,globalvalue*ratio]}. Has to be larger than 0 and smaller than 1. Defaults to 0.8.
#' @param seed seed for the random number generator used before creation of the landscape. Generator status will be saved and reset afterwards.
#'
#' @return returns a function. The function takes a point (vector) as input, with as many values as specified in \code{dimension}.
#' The function returns a vector value, representing the output from the two landscapes landscape.
#'
#' @references B. Yuan and M. Gallagher (2003) "On Building a Principled Framework for Evaluating and Testing Evolutionary Algorithms: A Continuous Landscape Generator". 
#' In Proceedings of the 2003 Congress on Evolutionary Computation, IEEE, pp. 451-458, Canberra, Australia.
#'
#' @seealso \code{\link{spotGlgCreate}} \code{\link{spotGlgCreateN}} \code{\link{spotGlgCreateRotSearched}}
#'
#' @export
###################################################################################
spotGlgCreateRot<- function(iDimension=2,alpha=pi/2,nGaussian=10,
				lower=rep(0,iDimension),upper=rep(1,iDimension),
				globalvalue=1,ratio=0.8,seed=1){
	fit <- spotGlgInit(iDimension,nGaussian,lower,upper,globalvalue,ratio,seed)
	fit$alpha = alpha
	fit$lower = lower
	fit$upper = upper
	fun <- function(x){
		spotGlgEvalRot(x,fit)$value
	}
	attributes(fun) <- fit
	fun
}

###################################################################################
#' Gaussian Landscape Evaluation (rotated)
#' 
#' This function is based on the Gaussian Landscape Generator by Bo Yuan and Marcus Gallagher.
#' It randomly evaluates one or several points in a combination of two landscapes, where the second is a rotation of the first.
#'
#' @param x matrix of sample sites, containing one point in each row.
#' @param glg list of values defining the Gaussian Landscape, created by \code{spotGlgInit}.
#'
#' @return returns a list, with the following items:\cr
#' \code{value} vector of 2 values, first belongs to non-rotated landscape.
#' \code{components} value of each component (will be NA)
#'
#' @references B. Yuan and M. Gallagher (2003) "On Building a Principled Framework for Evaluating and Testing Evolutionary Algorithms: A Continuous Landscape Generator". 
#' In Proceedings of the 2003 Congress on Evolutionary Computation, IEEE, pp. 451-458, Canberra, Australia.
#'
#' @seealso \code{\link{spotGlgCreateRot}}, \code{\link{spotGlgInit}}, \code{\link{spotGlgEval}}, \code{\link{spotGlgEvalN}}    
#'
#' @export
###################################################################################
spotGlgEvalRot <- function(x,glg){
	
	res <- spotGlgEval(x,glg)
	value <- res$value
	
	#mid of space:
	mid=((glg$upper-glg$lower)/2)+glg$lower
	x=x-mid	
	cosa=cos(glg$alpha)
	sina=sin(glg$alpha)
	rotation=cbind(c(cosa,sina),c(-sina,cosa))
	x=x%*%rotation
	x=x+mid
	
	
	res <- spotGlgEval(x,glg)
	value <- cbind(value, res$value)
	
	list(value=value,components=NA)
}


###################################################################################
#' Create Gaussian Landscape
#' 
#' This function is based on the Gaussian Landscape Generator by Bo Yuan and Marcus Gallagher.
#' It creates a new Gaussian Landscape every time it is called. This Landscape can be evaluated like a function.
#'
#' @param dimension dimensionality of the landscapes input space. Default is 2.
#' @param nGaussian number of Gaussian components in the landscape. Default is 10.
#' @param lower lower boundary of the landscape, defaults to \code{rep(0,dimension)}.
#' @param upper upper boundary of the landscape, defaults to \code{rep(1,dimension)}.
#' @param globalvalue the global maximum value, i.e. the maximum of the Gaussian component with the largest value. Default is 1.
#' @param ratio maximum ratio of the local maxima, local optima are randomly generated within \code{[0,globalvalue*ratio]}. Has to be larger than 0 and smaller than 1. Defaults to 0.8.
#' @param seed seed for the random number generator used before creation of the landscape. Generator status will be saved and reset afterwards.
#'
#' @return returns a function. The function takes a point (vector) as input, with as many values as specified in \code{dimension}.
#' The function returns a single scalar value, which is the Landscape value at the current point. The function has several attributes which reflect the items returned by \code{\link{spotGlgInit}}.
#' reflect the values returned by
#'
#' @references B. Yuan and M. Gallagher (2003) "On Building a Principled Framework for Evaluating and Testing Evolutionary Algorithms: A Continuous Landscape Generator". 
#' In Proceedings of the 2003 Congress on Evolutionary Computation, IEEE, pp. 451-458, Canberra, Australia.
#'
#' @examples 
#'	## Create a landscape function with default settings:
#'	landscapeFun <- spotGlgCreate()
#'	## Plot the landscape (uncomment before running this example)
#'	#spotSurf3d(landscapeFun)
#'
#' @seealso \code{\link{spotGlgCreateRot}} \code{\link{spotGlgCreateN}} \code{\link{spotGlgCreateRotSearched}}
#'
#' @export
###################################################################################
spotGlgCreate <- function(dimension=2,nGaussian=10,lower=rep(0,dimension),upper=rep(1,dimension),globalvalue=1,ratio=0.8,seed=1){
	fit <- spotGlgInit(dimension,nGaussian,lower,upper,globalvalue,ratio,seed)
	fun <- function(x){
		spotGlgEval(x,fit)$value
	}
	attributes(fun) <- fit
	fun
}


###################################################################################
#' Initialize Gaussian Landscape
#' 
#' This function is based on the Gaussian Landscape Generator by Bo Yuan and Marcus Gallagher.
#' It randomly initializes a Gaussian Landscape with the specified parameters.
#'
#' @param dimension dimensionality of the landscapes input space. Default is 2.
#' @param nGaussian number of Gaussian components in the landscape. Default is 10.
#' @param lower lower boundary of the landscape, defaults to \code{rep(0,dimension)}.
#' @param upper upper boundary of the landscape, defaults to \code{rep(1,dimension)}.
#' @param globalvalue the global maximum value, i.e. the maximum of the Gaussian component with the largest value. Default is 1.
#' @param ratio maximum ratio of the local maxima, local optima are randomly generated within \code{[0,globalvalue*ratio]}. Has to be larger than 0 and smaller than 1. Defaults to 0.8.
#' @param seed seed for the random number generator used before creation of the landscape. Generator status will be saved and reset afterwards.
#'
#' @return returns a list, with the following items:\cr
#' \code{mean} Matrix containing the mean vectors of the Gaussian components in the landscape, i.e. the locations of the local maxima of the functions. First vector (i.e. first row) will be the global maximum.
#' \code{covinv} Inverse of covariance matrix of each Gaussian component, stored as 3-dimensional array.
#' \code{opt} optimal values, i.e. maxima of the Gaussian components
#' \code{ngauss} number of Gaussian components
#' \code{d} is the \code{dimension}
#'
#' @references B. Yuan and M. Gallagher (2003) "On Building a Principled Framework for Evaluating and Testing Evolutionary Algorithms: A Continuous Landscape Generator". 
#' In Proceedings of the 2003 Congress on Evolutionary Computation, IEEE, pp. 451-458, Canberra, Australia.
#'
#' @author Original Matlab code by Bo Yuan, ported to R by Martin Zaefferer
#'
#' @seealso \code{\link{spotGlgCreate}}, \code{\link{spotGlgEval}} , \code{\link{spotGlgInitN}}
#'
#' @examples 
#'	## Create a landscape with default settings:
#'	landscape <- spotGlgInit()
#'	## Create a landscape with larger boundaries and more Gaussian components
#'	landscape <- spotGlgInit(2, 100, -5, 5, 10, 0.8)
#'
#' @export
###################################################################################
spotGlgInit<- function (dimension=2,nGaussian=10,lower=rep(0,dimension),upper=rep(1,dimension),globalvalue=1,ratio=0.8,seed=1){
	if(exists(as.character(substitute(.Random.seed))))
		SAVESEED<-.Random.seed
	else
		SAVESEED=NULL
	set.seed(seed)
	
	res <- spotGlgCore(dimension,nGaussian,lower,upper,globalvalue,ratio)

	if(!is.null(SAVESEED))
		assign(".Random.seed", SAVESEED, envir=globalenv())
	res
}

###################################################################################
#' Gaussian Landscape Core function
#' 
#' Core gaussian landscape function. Should not be called directly, as it does not contain proper seed handling.
#'
#' @param dimension dimensionality of the landscapes input space. Default is 2.
#' @param nGaussian number of Gaussian components in the landscape. Default is 10.
#' @param lower lower boundary of the landscape, defaults to \code{rep(0,dimension)}.
#' @param upper upper boundary of the landscape, defaults to \code{rep(1,dimension)}.
#' @param globalvalue the global maximum value, i.e. the maximum of the Gaussian component with the largest value. Default is 1.
#' @param ratio maximum ratio of the local maxima, local optima are randomly generated within \code{[0,globalvalue*ratio]}. Has to be larger than 0 and smaller than 1. Defaults to 0.8.
#'
#' @author Based on original Matlab code by Bo Yuan, ported to R and further adaptions by Martin Zaefferer
#'
#' @keywords internal
#' @export
###################################################################################
spotGlgCore <- function(dimension=2,nGaussian=10,lower=rep(0,dimension),upper=rep(1,dimension),globalvalue=1,ratio=0.8){
	if (dimension<=1|nGaussian<=0|any(upper<=lower)|globalvalue<=0|ratio<=0|ratio>=1){
		stop('Incorrect parameter values for gaussian landscape generator')}
	#% Generate rotation matrix

	e=diag(rep(1,dimension));   #% unit diagonal matrix

	rotation=array(0,c(dimension,dimension,nGaussian))
	covmatrix_inv <- rotation #initialize covmatrix
	
	for(i in 1:nGaussian){		
	   rotation[,,i]=e            #% initial rotation matrix for each Gaussian
	}

	for(i in 1:nGaussian){ 
		for(j in 1:(dimension-1)){        #% totally n(n-1)/2 rotation matrice
			for(k in (j+1):dimension){				  
				r=e
				
				alpha=runif(1)*pi/2-pi/4#% random rotation angle [-pi/4,pi/4]
				
				r[j,j]=cos(alpha)
				r[j,k]=sin(alpha)
				r[k,j]=-sin(alpha)
				r[k,k]=cos(alpha)

				rotation[,,i]=rotation[,,i]%*%r		
	}}}

	#% Generate covariance matrix

	variancerange=(upper-lower)/20  # % this controls the range of variances

	variance=matrix(runif(nGaussian*dimension),nGaussian,dimension)*variancerange+0.05;  #% add 0.05 to avoid zero variance     #TODO or mat mult?

	for(i in 1:nGaussian){
		  covmatrix=if(length(variance[i,])>1){diag(variance[i,])}else{variance[i,]}  #TODO similar in next line? 
		  covmatrix=if(length(variance[i,])>1){t(rotation[,,i])%*%covmatrix%*%rotation[,,i]}else{t(rotation[,,i])*covmatrix*rotation[,,i];}  #TODO similar in next line? 
		  covmatrix_inv[,,i]=solve(covmatrix)	  
	}

	#% Generate mean vectors randomly within [lower, upper]
	meanvector=matrix(runif(nGaussian*dimension),nGaussian,dimension)*(upper-lower)+lower   #TODO or mat mult?

	#% assign values to components
	optimumvalue=rep(0,nGaussian) #initialize
	optimumvalue[1]=globalvalue     #% the first Gaussian is set to be the global optimum

	#% values of others are randomly generated within [0,globalvalue*ratio]
	optimumvalue[2:nGaussian]=matrix(runif(1*(nGaussian-1)),1,nGaussian-1)*globalvalue*ratio
	list(mean=meanvector,covinv=covmatrix_inv,opt=optimumvalue,ngauss=nGaussian,d=dimension)
}

###################################################################################
#' Gaussian Landscape Evaluation
#' 
#' This function is based on the Gaussian Landscape Generator by Bo Yuan and Marcus Gallagher.
#' It randomly evaluates one or several points in a Gaussian Landscape created by \code{spotGlgInit}.
#'
#' @param x matrix of sample sites, containing one point in each row.
#' @param glg list of values defining the Gaussian Landscape, created by \code{spotGlgInit}.
#'
#' @return returns a list, with the following items:\cr
#' \code{value} value of the combined landscape
#' \code{components} value of each component
#'
#' @references B. Yuan and M. Gallagher (2003) "On Building a Principled Framework for Evaluating and Testing Evolutionary Algorithms: A Continuous Landscape Generator". 
#' In Proceedings of the 2003 Congress on Evolutionary Computation, IEEE, pp. 451-458, Canberra, Australia.
#'
#' @author Original Matlab code by Bo Yuan, ported to R by Martin Zaefferer
#'
#' @seealso \code{\link{spotGlgCreate}}, \code{\link{spotGlgInit}}, \code{\link{spotGlgEvalRot}}, \code{\link{spotGlgEvalN}} 
#'
#' @examples 
#'	## Create a landscape with default settings:
#'	landscape <- spotGlgInit()
#'	## Create a landscape with larger boundaries and more Gaussian components
#'	value <-  spotGlgEval(c(0.5,0.5),landscape)
#'
#' @export
###################################################################################
spotGlgEval <- function(x,glg){
	covmatrix_inv = glg$covinv #%the inverse covariance matrix of each component
	meanvector= glg$mean     #%the mean of each component
	optimumvalue= glg$opt   #%the peak value of each component
	nGaussian= glg$ngauss  #% total number of components
	n=glg$d #n: dimensionality
	
	p=nrow(x)                 #% p: number of individuals;      #TODO: this might be a problem, for one dimensional x ?
	if(is.null(p))p<-1	
	tmp=matrix(0,nGaussian,p)

	#%----------------------------------------------------
	for(i in 1:nGaussian) {             #% calculate the values generated by each component		
		newx=x-t(matrix(meanvector[i,],length(meanvector[i,]),p,byrow=FALSE))
		z=(newx%*%covmatrix_inv[,,i])*newx;			   
		tmp[i,]=rowSums(z);        	   
	}
	f=exp(-0.5*tmp/n);             #% f is a nGaussian-by-p matrix

	f=f*matrix(optimumvalue,length(optimumvalue),p,byrow=FALSE)#% multiply the peak value of each component
	                #% the value of each individual generated by each component
	value=apply(f,2,max)#max(f,[],1);      #% choose the maximum values as the fitness values
	list(value=value,components=f)
}

#glg<-initGLG (2, 100, -5, 5, 10, 0.8); 	#% 100 Gaussian components, [-5, 5]2, global optimum value 10, best possible local optimum value 8
#plotlandscapeGLG (-5, 5, 100,glg);  	#% Generate a surface plot and a contour plot with 100 samples within [-5, 5] in each dimension
#f=fitnessGLG(10*matrix(runif(2000),1000,2)-5,glg); 	#% Calculate the fitness values of 1000 random individuals within [-5, 5]2
