
setClass("BINOMresponse",contains="GLMresponse")

# method 'fit'
# use: in EM (M step)
# returns: (fitted) response with (new) estimates of parameters

# methods 'logDens' & dens
# use: instead of density slot in rModel
# returns: matrix with log(p(y|x,parameters))
setMethod("logDens","BINOMresponse",
	function(object) {
		if(NCOL(object@y) == 2) {
			dbinom(x=object@y[,1],size=rowSums(object@y),prob=predict(object),log=TRUE)
		} else {
			if(!NCOL(object@y==1)) stop("not a valid response matrix for BINOMresponse")
			dbinom(x=object@y,size=1,prob=predict(object),log=TRUE)
		}
	}
)

setMethod("dens","BINOMresponse",
	function(object,log=FALSE) {
		if(NCOL(object@y) == 2) {
			dbinom(x=object@y[,1],size=rowSums(object@y),prob=predict(object),log=log)
		} else {
			if(!NCOL(object@y==1)) stop("not a valid response matrix for BINOMresponse")
			dbinom(x=object@y,size=1,prob=predict(object),log=log)
		}
	}
)

setMethod("simulate",signature(object="BINOMresponse"),
	function(object,nsim=1,seed=NULL,times) {
    if(!is.null(seed)) set.seed(seed)
		if(missing(times)) {
			# draw in one go
			pr <- predict(object)
		} else {
			pr <- predict(object)[times,]
		}
		nt <- nrow(pr)
		n <- rowSums(object@y)
		if(NCOL(object@y) == 2) {
			response <- rbinom(nt*nsim,size=n,prob=pr)
 			response <- cbind(response,n-response)
		} else {
			response <- rbinom(nt*nsim,size=1,prob=pr)
		}
		#if(nsim > 1) response <- matrix(response,ncol=nsim)
 		response <- as.matrix(response)
		return(response)
	}
)
