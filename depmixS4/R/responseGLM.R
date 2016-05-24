
# 
# GLM response
# 

setClass("GLMresponse",
	representation(formula="formula",
		family="ANY"
	),
	prototype(
		formula=.~.,
		family=gaussian()
	),
	contains="response"
)


setMethod("GLMresponse",
	signature(formula="formula"),
	function(formula,data=NULL,family=gaussian(),pstart=NULL,fixed=NULL,prob=TRUE,na.action="na.pass", ...) {
		call <- match.call()
		mf <- match.call(expand.dots = FALSE)
		m <- match(c("formula", "data"), names(mf), 0)
		mf <- mf[c(1, m)]
		mf$drop.unused.levels <- TRUE
		mf$na.action <- na.action
		mf[[1]] <- as.name("model.frame")
		mf <- eval(mf, parent.frame())
		x <- model.matrix(attr(mf, "terms"),mf)
		if(any(is.na(x))) stop("'depmixS4' does not currently handle covariates with missing data.")
		y <- model.response(mf)
		if(!is.matrix(y)) y <- matrix(y,ncol=1)
		parameters <- list()
		constr <- NULL
		parameters$coefficients <- vector("numeric",length=ncol(x))
		names(parameters$coefficients) <- attr(x,"dimnames")[[2]]
		if(family$family=="gaussian") {
			parameters$sd <- 1
			constr <- list(
				parup = rep(Inf,ncol(x)+1),
				parlow = c(rep(-Inf,ncol(x)),0)
			)
		}
		if(family$family=="binomial") {
			# FIX ME
			y <- model.response(mf)
			if(NCOL(y) == 1) {
				if(is.factor(y)) y <- as.matrix(as.numeric(as.numeric(y)!=1)) else { # 21/06/12 changed this from "==" to "!=" in line with glm
					if(!is.numeric(y)) stop("model response not valid for binomial model")
					if(sum(y %in% c(0,1)) != length(y)) stop("model response not valid for binomial model")
					y <- as.matrix(y)
				}
			} else {
				if(ncol(y) != 2) {
					stop("model response not valid for binomial model")
				}
			}
		}
		if(family$family=="multinomial") {
			y <- model.response(mf)
			namesy <- NULL
			if(NCOL(y) == 1) {
				if(is.factor(y)) {
						namesy <- levels(y)
				    mf <- model.frame(~y-1,na.action=na.action)
				    y <- model.matrix(attr(mf, "terms"),mf)
					#y <- model.matrix(~y-1,na.action=na.action) 
				} else {
					if(!is.numeric(y)) stop("model response not valid for multinomial model")
					namesy <- levels(factor(y))
					mf <- model.frame(~factor(y)-1,na.action=na.action)
				    y <- model.matrix(attr(mf, "terms"),mf)
					#y <- model.matrix(~factor(y)-1,na.action=na.action)
				}
			}
			if(family$link=="mlogit") {
				if(any(y>1)) stop("multinomial response with n>1 not allowed with link='mlogit'")
				parameters$coefficients <- matrix(0,ncol=ncol(y),nrow=ncol(x))
				if(is.null(fixed)) {
					fixed <- parameters$coefficients
					fixed[,family$base] <- 1 
					fixed <- c(as.logical(t(fixed)))
			  }
			  if(is.null(namesy)) namesy <- 1:ncol(y)
			  colnames(parameters$coefficients) <- namesy
				rownames(parameters$coefficients) <- attr(x,"dimnames")[[2]]
# 				if(ncol(x)==1) names(parameters$coefficients) <- 1:ncol(y)
			}
			if(family$link=="identity") {
				if(ncol(x)>1) stop("covariates not allowed in multinomial model with identity link")
				ncy <- ncol(y)
				parameters$coefficients <- rep(1/ncy,ncy)
				if(is.null(namesy)) names(parameters$coefficients) <- paste("pr",1:ncol(y),sep="")
				else names(parameters$coefficients) <- namesy
				fixed <- rep(0,ncol(y)) 
				fixed <- c(as.logical(t(fixed)))
				constr <- list(
					lin = matrix(1,nrow=1,ncol=ncol(y)),
					linup = 1,
					linlow = 1,
					parup = rep(1,ncol(y)),
					parlow = rep(0,ncol(y))
			  )
		  }
		}
		npar <- length(unlist(parameters))
		if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
		if(!is.null(pstart)) {
				if(length(pstart)!=npar) stop("length of 'pstart' must be",npar)
				if(family$family=="multinomial") {
						if(family$link=="identity") {
								parameters$coefficients <- pstart[1:length(parameters$coefficients)]/sum(pstart[1:length(parameters$coefficients)])
						} else {
								if(prob) parameters$coefficients[1,] <- family$linkfun(pstart[1:ncol(parameters$coefficients)],base=family$base)
								else parameters$coefficients[1,] <- pstart[1:ncol(parameters$coefficients)]
						}
						pstart <- matrix(pstart,ncol(x),byrow=TRUE)
						if(ncol(x)>1) parameters$coefficients[2:ncol(x),] <- pstart[2:ncol(x),]
				} else {
						# if(prob) parameters$coefficients <- family$linkfun(as.numeric(pstart[1:length(parameters$coefficients)]))
						if(family$family=="binomial") {
								if(prob) parameters$coefficients[1] <- family$linkfun(pstart[1])
								else parameters$coefficients[1] <- pstart[1]
								if(ncol(x)>1) parameters$coefficients[2:ncol(x)] <- pstart[2:ncol(x)]
						} else {
								parameters$coefficients[1:ncol(x)] <- pstart[1:ncol(x)]
						}
						if(length(unlist(parameters))>length(parameters$coefficients)) {
								if(family$family=="gaussian") parameters$sd <- as.numeric(pstart[(length(parameters$coefficients)+1)])
						}
				}
		}
		mod <- switch(family$family,
				gaussian = new("NORMresponse",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar,constr=constr),
				binomial = new("BINOMresponse",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar,constr=constr),
				multinomial = new("MULTINOMresponse",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar,constr=constr),
				poisson = new("POISSONresponse",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar,constr=constr),
				Gamma = new("GAMMAresponse",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar,constr=constr),
				new("GLMresponse",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar,constr=constr)
		)
		mod
  }
)



setMethod("show","GLMresponse",
	function(object) {
		cat("Model of type ", object@family$family," (", object@family$link, "), formula: ", sep="")
		print(object@formula)
		cat("Coefficients: \n")
		print(object@parameters$coefficients)
		if(object@family$family=="multinomial"&object@family$link!="identity") {
			# also print probabilities at covariate values of zero
			cat("Probalities at zero values of the covariates.\n")
			if(!(is.null(dim(object@parameters$coefficients)))) {
				if(dim(object@parameters$coefficients)[1]>1) {
					cat(object@family$linkinv(object@parameters$coefficients[1,],base=object@family$base),"\n")
				} else {
					cat(object@family$linkinv(object@parameters$coefficients,base=object@family$base),"\n")
				}
			} else {
				if(object@family$link=="identity") cat(object@family$linkinv(object@parameters$coefficients),"\n")
				else {
					cat(object@family$linkinv(object@parameters$coefficients,base=object@family$base),"\n")
				}
			}
		}
		if(object@family$family=="binomial") {
			# also print probabilities at covariate values of zero
			cat("Probality at zero values of the covariates.","\n")
			cat(object@family$linkinv(object@parameters$coefficients[1]),"\n")
		}
		if(object@family$family=="gaussian") {
			cat("sd ",object@parameters$sd,"\n")
		}	
	}
)

setMethod("setpars","GLMresponse",
	function(object, values, which="pars", prob=FALSE, ...) {
		npar <- npar(object)
		if(length(values)!=npar) stop("length of 'values' must be",npar)
		# determine whether parameters or fixed constraints are being set
 		nms <- attributes(object@parameters$coefficients)
		if(length(values) == 0) return(object) # nothing to set; 
		switch(which,
			"pars"= {
				if(object@family$family=="multinomial") {
					object@parameters$coefficients <- matrix(values,ncol(object@x),byrow=TRUE)
					if(object@family$link=="mlogit") {	
						if(prob) object@parameters$coefficients[1,] <- object@family$linkfun(values[1:ncol(object@parameters$coefficients)],base=object@family$base)
					}
				} else {
					object@parameters$coefficients <- values[1:length(object@parameters$coefficients)] # matrix(values,ncol(object@x),byrow=TRUE) # this needs fixing!!!!
				}
				if(length(unlist(object@parameters))>length(object@parameters$coefficients)) {
					if(object@family$family=="gaussian") object@parameters$sd <- values[(length(object@parameters$coefficients)+1)]
				}
			},
			"fixed" = {
				object@fixed <- as.logical(values)
			}
	  )
	  attributes(object@parameters$coefficients) <- nms
		return(object)
	}
)

setMethod("getpars","GLMresponse",
	function(object,which="pars",...) {
		switch(which,
			"pars" = {
				parameters <- numeric()
				if(object@family$family=="multinomial"&object@family$link=="mlogit") {
					# coefficient is usually a matrix here
					tmp <- object@parameters$coefficients
					parameters <- c(t(tmp)) # Why transpose?
					names(parameters) <- paste(rep(rownames(tmp),each=length(colnames(tmp))),colnames(tmp),sep=".")
				} else {
					parameters <- object@parameters$coefficients
					if(object@family$family=="gaussian") {
						nms <- names(parameters)
						parameters <- c(parameters,object@parameters$sd)
						names(parameters) <- c(nms,"sd")
					}
					
				}
				pars <- parameters
			},
			"fixed" = {
				pars <- object@fixed
			}
		)
		return(pars)
	}
)

# methods: fit, logDens, predict
# use: in EM (M step)
# returns: (fitted) response with (new) estimates of parameters

setMethod("fit","GLMresponse",
	function(object,w) {
    if(missing(w)) w <- NULL
		pars <- object@parameters
		start <- pars$coefficients
		start[is.na(start)] <- 0
		fit <- glm.fit(x=object@x,y=object@y,weights=w,family=object@family,start=start)
		pars$coefficients <- fit$coefficients
		object <- setpars(object,unlist(pars))
		object
	}
)

setMethod("logLik","GLMresponse",
	function(object) {
		sum(logDens(object))
	}
)

setMethod("predict","GLMresponse",
	function(object) {
	  nas <- is.na(object@parameters$coefficients)
    if(sum(nas) == 0) {
      object@family$linkinv(object@x%*%object@parameters$coefficients)
    } else {
		  object@family$linkinv(as.matrix(object@x[,!nas])%*%object@parameters$coefficients[!nas])
    }
	}
)
