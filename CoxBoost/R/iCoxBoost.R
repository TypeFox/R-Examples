cvcb.control <- function(K=10,type=c("verweij","naive"),parallel=FALSE,upload.x=TRUE,multicore=FALSE,
                 folds=NULL) 
{
	type <- match.arg(type)
	list(K=K,type=type,parallel=parallel,upload.x=upload.x,multicore=multicore,folds=folds)
}

iCoxBoost <- function(formula,data=NULL,weights=NULL,subset=NULL,
					  mandatory=NULL,cause=1,standardize=TRUE,
					  stepno=200,criterion=c("pscore","score","hpscore","hscore"),
					  nu=0.1,stepsize.factor=1,varlink=NULL,
					  cv=cvcb.control(),trace=FALSE,...)
{
	criterion <- match.arg(criterion)

    call <- match.call(expand.dots=FALSE)
    formula <- eval(call$formula)
    
    if (is.null(data)) {
    	actual.terms <- terms(formula)
    	response <- model.response(model.frame(formula,data=environment(formula)))
	    x <- model.matrix(actual.terms)[,-c(1),drop=FALSE]
    } else {
    	actual.terms <- terms(formula,data=data)
    	response <- model.response(model.frame(formula,data))
	    x <- model.matrix(actual.terms,data=data)[,-c(1),drop=FALSE]
    }

    actual.time <- as.numeric(response[,"time"])
    if (class(response) == "Hist") {
        actual.status <- rep(2,NROW(response))
        actual.status[as.numeric(response[,"event"] == cause) == 1] <- 1
        actual.status[as.numeric(response[,"status"]) == 0] <- 0
    } else {
    	actual.status <- as.numeric(response[,"status"])    	
    }

    if (is.null(subset)) subset <- 1:length(actual.time)
    penalty <- sum(actual.status[subset]==cause)*(1/nu-1)

    actual.names <- colnames(x)
    unpen.index <- NULL
    if (!is.null(mandatory)) {
    	if (sum(mandatory %in% actual.names) < length(mandatory)) {
    		stop("Some variables in 'mandatory' are not part of the model.")
    	}
    	unpen.index <- which(actual.names %in% mandatory)
    }

    pendistmat <- NULL
    connected.index <- NULL
    if (!is.null(varlink)) {
    	if (!is.list(varlink) || length(varlink) < 2) {
    		stop("'varlink' has to be a list with at least two elements")
    	}

    	actual.source <- varlink[[1]]
    	actual.target <- varlink[[2]]

    	if (length(varlink) == 2) {
    		actual.factor <- rep(1,length(actual.source))
    	} else {
			actual.factor <- varlink[[3]]    		
    	}

    	if (length(actual.source) != length(actual.target) || length(actual.source) != length(actual.factor)) {
    		stop("Source and target (and factor) vectors of 'varlink' have to be of the same length.")
    	}

    	actual.connected <- union(actual.source,actual.target)
    	if (!all(actual.connected %in% actual.names)) {
    		stop("Some elements of 'varlink' do not appear to be part of the model.")
    	}

    	connected.index <- which(actual.names %in% actual.connected)
    	actual.source.match <- match(actual.source,actual.names[connected.index])
    	actual.target.match <- match(actual.target,actual.names[connected.index])

    	pendistmat <- Matrix(0,length(connected.index),length(connected.index))

    	for (i in seq(along=actual.source.match)) {
    		pendistmat[actual.source.match[i],actual.target.match[i]] <- actual.factor[i]
    	}
    }

    use.stepno <- stepno
    cv.res <- NULL
    if (is.list(cv) || cv == TRUE) {
    	if (!is.list(cv)) cv <- cvcb.control()
	    cv.res <- cv.CoxBoost(time=actual.time,status=actual.status,x=x,
	    					  maxstepno=stepno,K=cv$K,type=cv$type,parallel=cv$parallel,
	    					  upload.x=cv$upload.x,multicore=cv$multicore,folds=cv$folds,
	    					  weights=weights,subset=subset,standardize=standardize,
	    					  penalty=penalty,criterion=criterion,
	    					  stepsize.factor=stepsize.factor,
	    					  unpen.index=unpen.index,pendistmat=pendistmat,
	    					  connected.index=connected.index,
	    					  trace=trace,...)
	    use.stepno <- cv.res$optimal.step
    }

    res <- CoxBoost(time=actual.time,status=actual.status,x=x,
    				weights=weights,subset=subset,standardize=standardize,
    				penalty=penalty,criterion=criterion,
    				stepno=use.stepno,unpen.index=unpen.index,pendistmat=pendistmat,
	    			stepsize.factor=stepsize.factor,
    				connected.index=connected.index,
    				trace=trace,...)
 
    res$call <- call
    res$formula <- formula
    res$terms <- actual.terms
    res$cause <- cause
    res$cv.res <- cv.res
    
    class(res) <- c("iCoxBoost",class(res))
    res
}

predict.iCoxBoost <- function(object,newdata=NULL,subset=NULL,at.step=NULL,times=NULL,
							  type=c("lp","logplik","risk","CIF"),...) 
{
    type <- match.arg(type)

	newtime <- NULL
	newstatus <- NULL

	if (is.null(newdata)) {
		newdata <- NULL
	} else {
		if (type == "logplik") {
			new.response <- model.response(model.frame(object$formula,newdata))
			newtime <- as.numeric(new.response[,"time"])

			if (class(new.response) == "Hist") {
                newstatus <- rep(2,NROW(new.response))
                newstatus[as.numeric(new.response[,"event"] == object$cause) == 1] <- 1
                newstatus[as.numeric(new.response[,"status"]) == 0] <- 0
			} else {
				newstatus <- as.numeric(new.response[,"status"])    	
			}
		}

		newdata <- model.matrix(object$terms,
				   data=model.frame(object$formula,data=newdata))[,-c(1),drop=FALSE]
	}

	predict.CoxBoost(object,newdata,newtime,newstatus,subset,at.step,times,type,...)
}

