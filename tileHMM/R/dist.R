# dist.R
# 
# S4 classes and methods representing distribution objects
#
# Author: Peter Humburg
###############################################################################

## A virtual class representing a distribution
setClass("dist",representation(),prototype())

## New generic functions
setGeneric("sampleObs",def=function(dist,size,...){standardGeneric("sampleObs")})

## A class representing a discrete distribution over a given alphabet
setClass("discDist",representation(alpha="character",prob="numeric",dstr="list"),
		prototype(alpha=character(0),prob=numeric(0),dstr=list()), contains="dist")

## Initialising discDist
## Expects either two vectors specifying alphabet and associated probabilities
## or a list with named entries. Names correspond to the alphabet, entries to probabilities.
setMethod("initialize","discDist",
		function(.Object,alpha=character(0),prob=numeric(0),dstr=list()){
			if(length(dstr) > 0){
				alpha <- names(dstr)
				prob <- as.numeric(dstr)
			}
			
			if(length(alpha) != length(prob)) stop("Specified alphabet and distribution have different length!")
			psum <- sum(prob)
			pmin <- min(prob)
			if(pmin < 0) stop(paste("Probabilities have to be greater than 0, found ",pmin))
			if(psum == 0) stop("Illegal probability distribution: sum of probabilities is 0!")
			if(psum + 1 != 2){
				warning(paste("Illegal probability distribution: sum of probabilities is ",psum,".",sep=''))
				prob <- prob/psum
				warning("Probabilities were normalised.")
			}
			if(length(dstr) == 0){
				dstr <- list()
				for(i in 1:length(alpha)) dstr[[alpha[i]]] <- prob[i] 
			}
			
			.Object@alpha <- alpha
			.Object@prob <- prob
			.Object@dstr <- dstr
			.Object
		}
)
## accessing a discDist object
"[.discDist" <- function(x,i,...){x@dstr[[i]]}
"[[.discDist" <- function(x,i,...){x@dstr[[i]]}


## converting discDist objects to simple data structures
## data.frame
as.data.frame.discDist <- function(x,row.names=NULL,optional=FALSE,...){
	as.data.frame(x@dstr,row.names,optional)
}

## matrix
as.matrix.discDist <- function(x,...){
	t(as.matrix(x@dstr))
}
	
## vector
as.vector.discDist <- function(x,mode="any"){
	vec <- as.vector(x@prob,mode)
	names(vec) <- x@alpha
	vec
}

## define length of a discrete distribution as the number of elements
setMethod("length","discDist",function(x) length(x@alpha))

## Get a sample from the distribution
setMethod("sampleObs",c("discDist","numeric"),
		function(dist,size){
			obs <- sample(dist@alpha,size,T,prob=dist@prob)
			obs
		}
)

## A super class for continuous distributions
setClass("contDist",representation(components="matrix"),
		prototype(components=matrix(ncol=3,nrow=0)), contains="dist")
setMethod("initialize","contDist",
	function(.Object,weight=1,center=0,disp=1){
		## check parameters for consistency
		length.w <- length(weight)
		length.c <- length(center)
		length.d <- length(disp)
		if(length.w != length.c || length.w != length.d || length.c != length.d) {
			stop("length of parameters has to agree!")
		}
		if(min(weight) < 0){
			stop("Mixture weights have to be positive!")
		}
		if(sum(weight) != 1){
			warning("Sum of mixture weights is not 1. Weights will be normalised.")
			weight <- weight/sum(weight)
		}
		components <- cbind(weight,center,disp)
		colnames(components) <- c("weight","centrality","dispersion")
		.Object@components <- components
		.Object
	}
)

## t distribution
setClass("tDist", representation(), prototype(), contains="contDist")
setMethod("initialize","tDist",
	function(.Object,mean=0,var=1,df=3){
		.Object <- callNextMethod(.Object,weight=1,center=mean,disp=var)
		.Object@components <- cbind(.Object@components,df)
		colnames(.Object@components) <- c("weight","mean","variance","df")
		.Object
	}
)

#############################################################################
## Methods to access continuous density distributions and derived classes  ##
#############################################################################

## returns selected parameter of i-th mixture component
## selected parameter of all mixture components if j is missing
"[.contDist" <- function(x,i,j,...){
	if(missing(i)) value <- x@components[,j] else value <- x@components[i,j]
	value
}

## changing mixture components
"[<-.contDist" <- function(x,i,j,value){
	if(missing(i) && !missing(j)){
		n <- dim(x@components)[1]
		x[1:n,j] <- value
	} 
	if(!missing(i) && missing(j)){
		n <- dim(x@components)[2]
		x[i,1:n] <- value
	}
	if(missing(i) && missing(j)){
		n <- dim(x@components)[1]
		m <- dim(x@components)[2]
		x[1:n,1:m] <- value
	}
	if(!missing(i) && !missing(j)){
		x@components[i,j] <- value
		if(class(x@components[,c(1,3)]) == "matrix"){
			m <- apply(x@components[,c(1,3)],2,min)
		} else{
			m <- c(x@components[,c(1,3)])
		}		
		if(m[1] < 0){
			stop("Weights have to be positive!")
		}
		if(m[2] < 0){
			stop("Variance has to be positive!")
		}
		if(sum(x@components[,1]) != 1){
			warning("Weights will be renormalised.")
			x@components[,1] <- x@components[,1]/sum(x@components[,1])
		}
	}
	x
}

"[<-.tDist" <- function(x,i,j,value){"[<-.contDist"(x,i,j,value)}

## tDist specific methods
## returns weighted density of i-th mixture component at point j if j is numeric
## or weighted density of all mixture components at point j if j is numeric and i is missing
## If j is a character string the corresponding parameter value is returned
"[.tDist" <- function(x,i,j,log=FALSE,...){
	if(!missing(i)){
		if(is(j,"numeric")){
			ret <- dt((j-x@components[i,2])/sqrt(x@components[i,3]),df=x@components[i,4],log=log)
			if(log){
				ret <- log(x@components[i,1]) + ret - 0.5 * log(x@components[i,3])
			} else {
				ret <- x@components[i,1] * ret / sqrt(x@components[i,3])
			}
		} else{
			ret <- "[.contDist"(x,i,j,...)
		}
	} else{
		if(is(j,"numeric")){
			ret <- dt((j-x@components[,2])/sqrt(x@components[,3]),df=x@components[,4],log=log)
			if(log){
				ret <- log(x@components[,1]) + ret - 0.5 * log(x@components[,3])
			} else{
				ret <- x@components[,1] * ret / sqrt(x@components[,3])
			}
		} else{
			ret <- "[.contDist"(x, ,j,...)
		}
	}
	ret
}


#######################################################
## Plotting and printing contDist objects            ##
#######################################################
plot.contDist <- function(x, step.size=0.01, new.plot=TRUE, weight=1, ...){
	xrange <- list(...)$xlim
		
	if(is.null(xrange)){
		min.mean <- which.min(x@components[,2])
		max.mean <- which.max(x@components[,2])
		xrange <- c(x@components[min.mean,2]-4*sqrt(x@components[min.mean,3]),
					x@components[max.mean,2]+4*sqrt(x@components[max.mean,3]))
	}
	
	sample <- seq(xrange[1],xrange[2],by=step.size)
	dcomp <- function(obs,d){
		dens <- d[,obs]
		dens
	}
	dmix <- as.matrix(apply(as.matrix(sample),1,dcomp,x))
	if(dim(x@components)[1] > 1) dmix <- t(dmix)
	dmix <- apply(dmix,1,sum) * weight
	if(new.plot){
		plot(sample,dmix,type="l",xlab="",ylab="Density",...)
		abline(h=0)
	} else{
		lines(sample,dmix,...)
	}
	invisible(dmix)
}

plot.tDist <- plot.contDist

setMethod("show", "contDist",
	function(object){
		cat("An object of class \"", class(object), "\"\n", sep='')
		print(object@components)
	}
)



#########################################################
## Sample from contDist objects                        ##
#########################################################
## Get a sample from a t distribution
setMethod("sampleObs",signature=c("tDist","numeric"),
	function(dist,size){
		n <- dim(dist@components)[1]
		comp <- sample(1:n,size,T,prob=dist@components[,1])
		obs <- numeric()
		for(i in comp){
			obs <- c(obs,(rt(1,df=dist@components[i,4])+dist@components[i,2])*sqrt(dist@components[i,3]))
		}
		obs
	}
)


