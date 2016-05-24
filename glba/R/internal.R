pmod <-
function(formula, data, prefix) {
	mf <- model.frame(formula=formula,data=data)
	x <- model.matrix(formula,data=mf)
	pars <- numeric(ncol(x))
	names(pars) <- paste(prefix,dimnames(x)[[2]])
	ret <- list(x=x,pars=pars,formula=formula,prefix=prefix)
	class(ret) = "pmod"
	return(ret)
}



getPars <-
function(object, ...) {
	return(object$pars)
}

setPars <-
function(object,values, ...) {
	np <- names(object$pars)
	object$pars <- values
	names(object$pars) <- np
	return(object)
}

summary.pmod <-
function(object, ...) {
  	cat("\nmodel (", object$prefix,"): ",sep="")
	print(object$formula, showEnv=FALSE)
 	cat("parameters\n")
 	cat(object$pars, "\n")
}

print.pmod <-
function(x, ...) {
  	cat("\nmodel (", x$prefix,"): ",sep="")
	print(x$formula, showEnv=FALSE)
 	cat("parameters\n")
 	cat(x$pars, "\n")
}

predpmod <-
function(object,pars) {return(object$x%*%pars)}

predict.pmod <-
function(object, ...) {return(object$x%*%object$pars)}

reorderDrift <-
function(resp,drift) {
	nresp <- ncol(drift)
	nn <- nrow(drift)
 	if(all.equal(sort(unique(resp)),c(0,1))) resp <- ifelse(resp==0,2,1)
 	else resp <- as.numeric(as.factor(resp))
	resp <- as.numeric(as.factor(resp))
	ordResp <- function(x, nr){ c(x,(1:nr)[-x])}
	ord <- t(sapply(resp,ordResp,nr=nresp))
	dot <- matrix(,ncol=nresp,nrow=nn)
	for(i in 1:ncol(ord)) {
		ind <- matrix(c(1:nn,ord[,i]),ncol=2)
		dot[,i] <- drift[ind]
	}
	dot
}

obj <-
function(rt,pars,loglink,weights) {
	# vectorized loglike function
	# rt: a vector with response times
	# pars: matrix with 4+nrcat parameters on each row to model each rt
	# the drift pars are ordered: the drift for the given response first, the others 
	# after that (order in the remaining drifts does not make a difference)
	for(i in 1:4) if(loglink[i]) pars[,i]=exp(pars[,i])
	ndrift <- dim(pars)[2]-4
	if(ndrift<2) stop("nr of drift pars should at least be two")
	ll <- numeric(length(rt))
	
	ll <- n1PDF(t=rt-pars[,4], x0max=pars[,2],
		chi=pars[,2]+pars[,3], sdI=pars[,1], # sdI=0.15, # Scott: I fit chi-x0max.
		drift=pars[,5:(4+ndrift)])	
	
# 	return(logl=-sum(log(pmax(weights*ll,1e-10)))) # this has weird effects due to the contaminant model ...
	return(logl=sum(log(weights*ll)))
}
