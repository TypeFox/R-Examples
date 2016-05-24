# 
# for the transition models and the prior (y is missing, ie there is no
# response, and nstates must be provided as the number of categories
# neccessary in the mulinomial model)
# 

setClass("transInit",contains="GLMresponse")

setMethod("transInit",
	signature(formula="formula"),
	function(formula,nstates,data=NULL,family=multinomial(),pstart=NULL,fixed=NULL,prob=TRUE, ...) {
		call <- match.call()
		if(formula==formula(~1) &is.null(data)) { # &is.null(data) removed this in the condition as it 
			# creates the wrong dimension for the dens function when data is used (but not needed)
			x <- matrix(1,ncol=1)
		} else {
			mf <- match.call(expand.dots = FALSE)
			m <- match(c("formula", "data"), names(mf), 0)
			mf <- mf[c(1, m)]
			mf$drop.unused.levels <- TRUE
			mf[[1]] <- as.name("model.frame")
			mf <- eval(mf, parent.frame())
			x <- model.matrix(attr(mf, "terms"),mf)
			if(any(is.na(x))) stop("'depmixS4' does not currently handle covariates with missing data.")
		}
		y <- matrix(1,ncol=1) # y is not needed in the transition and init models
		parameters <- list()
		constr <- NULL
		if(is.null(nstates)) stop("'nstates' must be provided in call to transInit model")
		if(family$family=="multinomial") {
			if(family$link=="identity") {
					if(ncol(x)>1) stop("covariates not allowed in multinomial model with identity link")
					parameters$coefficients <- rep(1/nstates,nstates)
					names(parameters$coefficients) <- paste("pr",1:nstates,sep="")
				if(is.null(fixed)) {
					fixed <- matrix(0,nrow=1,ncol=nstates)
					fixed <- rep(0,nstates) # this needs to be fixed at some point using contraints
					fixed <- c(as.logical(fixed))
				}
				constr <- list(
					lin = matrix(1,nrow=1,ncol=nstates),
					linup = 1,
					linlow = 1,
					parup = rep(1,nstates),
					parlow = rep(0,nstates)
				)
			} else {
				parameters$coefficients <- matrix(0,ncol=nstates,nrow=ncol(x))
				if(is.null(fixed)) {
					fixed <- parameters$coefficients
					fixed[,family$base] <- 1 
					fixed <- c(as.logical(t(fixed)))
			  }
				colnames(parameters$coefficients) <- paste("St",1:nstates,sep="")
				rownames(parameters$coefficients) <- attr(x,"dimnames")[[2]]
			}
		}
		npar <- length(unlist(parameters))
		if(is.null(fixed)) fixed <- rep(0,npar)
		if(!is.null(pstart)) {
			if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar)
			if(family$family=="multinomial") {
				if(family$link=="identity") {
					parameters$coefficients[1:nstates] <- pstart[1:nstates]
					parameters$coefficients <- parameters$coefficients/sum(parameters$coefficients)
				} else {
					if(prob) {
						parameters$coefficients[1,] <- family$linkfun(pstart[1:ncol(parameters$coefficients)],base=family$base)
					} else {
						parameters$coefficients[1,] <- pstart[1:ncol(parameters$coefficients)]
					}
					pstart <- matrix(pstart,ncol(x),byrow=TRUE)
					if(ncol(x)>1) parameters$coefficients[2:ncol(x),] <- pstart[2:ncol(x),]
				}
			} else {
				if(family$link=="identity") parameters$coefficients <- family$linkfun(pstart[1:length(parameters$coefficients)])
				else parameters$coefficients <- family$linkfun(pstart[1:length(parameters$coefficients)],base=family$base)
			}
		}
		# FIX this: do we need a switch here?
		mod <- switch(family$family,
			multinomial = new("transInit",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar,constr=constr),
			new("transInit",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar,constr=constr)
		)
		mod
	}
)

setMethod("logDens","transInit",
	function(object) {
		log(predict(object))
	}
)

setMethod("dens","transInit",
	function(object,log=FALSE) {
		if(log) log(predict(object))
		else predict(object)
	}
)

setMethod("predict","transInit",
	function(object,newx) {
		if(missing(newx)) {
			if(object@family$link=="identity") object@family$linkinv(object@x%*%object@parameters$coefficients)
			else object@family$linkinv(object@x%*%object@parameters$coefficients,base=object@family$base)
		} else {
			if(!(is.matrix(newx))) stop("'newx' must be matrix in predict(transInit)")
			if(!(ncol(newx)==nrow(object@parameters$coefficients))) stop("Incorrect dimension of 'newx' in predict(transInit)")
			if(object@family$link=="identity") object@family$linkinv(newx%*%object@parameters$coefficients)
			else object@family$linkinv(newx%*%object@parameters$coefficients,base=object@family$base)
		}
	}
)

setMethod("fit","transInit",
	function(object,w,ntimes) {
		pars <- object@parameters
		if(missing(w)) w <- NULL
		pars <- object@parameters
		base <- object@family$base # delete me
		#y <- object@y[,-base]
		y <- object@y
		x <- object@x
		if(is.matrix(y)) na <- unlist(apply(y,2,function(x) which(is.na(x)))) else na <- which(is.na(y))
		if(is.matrix(x)) na <- c(na,unlist(apply(x,2,function(x) which(is.na(x))))) else na <- c(na,which(is.na(x)))
		if(!is.null(w)) na <- c(na,which(is.na(w)))
		y <- as.matrix(y)
		x <- as.matrix(x)
		na <- unique(na)
		if(length(na)>0) {
			x <- x[-na,]
			y <- y[-na,]
			#y <- round(y) # delete me
			if(!is.null(w)) w <- w[-na]
		}
		switch(object@family$link,
		  mlogit = {
    		mask <- matrix(1,nrow=nrow(pars$coefficients),ncol=ncol(pars$coefficients))
    		mask[,base] <- 0 # fix base category coefficients to 0
    		mask <- rbind(0,mask) # fix "bias" nodes to 0
    		Wts <- mask
    		Wts[-1,] <- pars$coefficients # set starting weights
    		Wts[Wts == Inf] <- .Machine$double.max.exp # Fix this!!!!
    		Wts[Wts == -Inf] <- .Machine$double.min.exp # Fix this!!!!!
    		if(!is.null(w)) {
    			if(NCOL(y) < 3) {
    				fit <- nnet.default(x,y,weights=w,size=0,entropy=TRUE,skip=TRUE,mask=mask,Wts=Wts,trace=FALSE)
    			} else {
    				fit <- nnet.default(x,y,weights=w,size=0,softmax=TRUE,skip=TRUE,mask=mask,Wts=Wts,trace=FALSE)
    			}
    		} else {
    			if(NCOL(y) < 3) {
    				fit <- nnet.default(x,y,size=0,entropy=TRUE,skip=TRUE,mask=mask,Wts=Wts,trace=FALSE)
    			} else {
    				fit <- nnet.default(x,y,size=0,softmax=TRUE,skip=TRUE,mask=mask,Wts=Wts,trace=FALSE)
    			}
    		}
    		pars$coefficients <- t(matrix(fit$wts,ncol=ncol(pars$coefficients),nrow=nrow(pars$coefficients)+1)[-1,])
    		object <- setpars(object,unlist(pars))
		},
		identity = {
				if(!is.null(w)) {
						sw <- sum(w)
						pars <- colSums(w*object@y)/sum(w)
				} else {
						pars <- colMeans(object@y)
				}
 				pars[pars<1e-6] <- 0 # set small values to zero
				pars <- pars/sum(pars)
				object <- setpars(object,pars)
		},
		stop("link function not implemented")
	)
	object
}
)

setMethod("simulate",signature(object="transInit"),
	function(object,nsim=1,seed=NULL,times,is.prior=FALSE,...) {
		if(!is.null(seed)) set.seed(seed)
		if(is.prior) {
			pr <- dens(object)
			sims <- array(apply(pr,1,rmultinom,n=nsim,size=1),dim=c(ncol(pr),nsim,nrow(pr)))
			states <- t(apply(sims,c(2,3), function(x) which(x==1)))
			return(states)
		} else {
			if(missing(times)) {
				# this is likely to be a homogeneous model...???
				pr <- predict(object)
			} else {
				pr <- predict(object)[times,]
				if(length(times)==1) pr <- matrix(pr,ncol=length(pr))
			}
			nt <- nrow(pr)
			sims <- array(apply(pr,1,rmultinom,n=nsim,size=1),dim=c(ncol(pr),nsim,nt))
			states <- t(apply(sims,c(2,3), function(x) which(x==1)))
			# states <- apply(apply(pr,2,rmultinom rmultinom(nt*nsim,size=1,prob=pr),2,function(x) which(x==1))
			return(states)
		}
	}
)

setMethod("getdf","transInit",
	function(object) {
		npar <- sum(!object@fixed)
		if(object@family$link == "identity") {
			npar <- npar-1
			if(npar<0) npar <- 0
		}
		npar
	}
)
