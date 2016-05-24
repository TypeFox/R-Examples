
setClass("MULTINOMresponse",contains="GLMresponse")

setMethod("fit","MULTINOMresponse",
	function(object,w) {
		if(missing(w)) w <- NULL
		nas <- is.na(rowSums(object@y))
		if(object@family$link=="mlogit") {
			pars <- object@parameters
			base <- object@family$base # delete me
			y <- as.matrix(object@y[!nas,])
			x <- as.matrix(object@x[!nas,])
			#if(is.null(w)) w <- rep(1,nrow(y))
			# mask is an nx*ny matrix (x are inputs, y are output levels)
			mask <- matrix(1,nrow=nrow(pars$coefficients),ncol=ncol(pars$coefficients))
			mask[,base] <- 0 # fix base category coefficients to 0
			mask <- rbind(0,mask) # fix automatic "bias" nodes to 0
			Wts <- mask
			Wts[-1,] <- pars$coefficients # set starting weights
			if(!is.null(w)) {
				if(NCOL(y) < 3) {
					fit <- nnet.default(x,y,weights=w[!nas],size=0,entropy=TRUE,skip=TRUE,mask=mask,Wts=Wts,trace=FALSE)
				} else {
					fit <- nnet.default(x,y,weights=w[!nas],size=0,softmax=TRUE,skip=TRUE,mask=mask,Wts=Wts,trace=FALSE)
				}
			} else {
				if(NCOL(y) < 3) {
					fit <- nnet.default(x,y,size=0,entropy=TRUE,skip=TRUE,mask=mask,Wts=Wts,trace=FALSE)
				} else {
					fit <- nnet.default(x,y,size=0,softmax=TRUE,skip=TRUE,mask=mask,Wts=Wts,trace=FALSE)
				}
			}
			# this is necessary because setpars wants coefficients in column major order
			pars$coefficients <- t(matrix(fit$wts,ncol=ncol(pars$coefficients),nrow=nrow(pars$coefficients)+1)[-1,])
			# parameters are set correctly now
			object <- setpars(object,unlist(pars))
		}
		if(object@family$link=="identity") {
				if(is.null(w)) w <- rep(1,nrow(object@y))
				sw <- sum(w[!nas])
				pars <- c(apply(as.matrix(object@y[!nas,]),2,function(x){sum(x*w[!nas])/sw}))
# 				if(any(pars<1e-5)) warning("Parameters smaller than 1e-5 have been set to zero.")
				if(!all(pars < 1e-6)) pars[pars<1e-6] <- 0 # set small values to zero
				pars <- pars/sum(pars)
				object <- setpars(object,pars)
		}
		object
	}
)

setMethod("logDens","MULTINOMresponse",
	function(object) {
	    rsums <- rowSums(object@y)
		if(all(rsums[!is.na(rsums)]==1)) {
			return(log(rowSums(object@y*predict(object))))
		} else {
			nr <- nrow(object@y)
			res <- matrix(nrow=nr)
			pr <- predict(object)
			# fix this loop!!!! replace with call to apply? or dmultinomial? or other vectorized version?
			# possibly use dmultinomial in package mcd2
			for(i in 1:nrow(object@y)) {
				res[i,1] <- dmultinom(object@y[i,],prob=pr[i,])
			}
			return(log(res))
		}
	}
)

setMethod("dens","MULTINOMresponse",
	function(object,log=FALSE) {
	    rsums <- rowSums(object@y)
		if(all(rsums[!is.na(rsums)]==1)) {
			if(log) return(log(rowSums(object@y*predict(object))))
			else return(rowSums(object@y*predict(object)))
		} else {
			nr <- nrow(object@y)
			res <- matrix(nrow=nr)
			pr <- predict(object)
			# fix this loop!!!! replace with call to apply? or dmultinomial? or other vectorized version?
			# possibly use dmultinomial in package mcd2
			for(i in 1:nrow(object@y)) {
				res[i,1] <- dmultinom(object@y[i,],prob=pr[i,])
			}
			if(log) return(log(res)) 
			else return(res)
		}
	}
)

setMethod("predict","MULTINOMresponse",
	function(object) {
		if(object@family$link=="identity") object@x%*%object@parameters$coefficients
		else {
			object@family$linkinv(object@x%*%object@parameters$coefficients,base=object@family$base)
		}
	}
)

setMethod("simulate",signature(object="MULTINOMresponse"),
	function(object,nsim=1,seed=NULL,times) {
		if(!is.null(seed)) set.seed(seed)
		if(missing(times)) {
			# draw all times in one go
			pr <- predict(object)
			n <- rowSums(object@y)
		} else {
			pr <- predict(object)[times,]
			n <- rowSums(object@y)[times]
			if(length(times)==1) pr <- matrix(pr,ncol=length(pr))
		}
		nt <- nrow(pr)
		sims <- array(apply(pr,1,rmultinom,n=nsim,size=n),dim=c(ncol(pr),nsim,nt))
		sims <- matrix(aperm(sims,c(3,2,1)),nrow=nsim*nt,ncol=ncol(pr))
		#response <- t(apply(sims,c(2,3), function(x) which(x==1)))
		return(sims)
	}
)

setMethod("getdf","MULTINOMresponse",
	function(object) {
		npar <- sum(!object@fixed)
		if(object@family$link == "identity") {
			npar <- npar-1
			if(npar<0) npar <- 0
		}
		npar
	}
)
