metatest <-
function(formula, variance, data, threshold=0.00001, maxiter=100, npermut=1000, ...) {

	call <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula","variance", "data"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	x <- model.matrix(attr(mf, "terms"),mf)
	y <- model.response(mf)
	if(!is.matrix(y)) y <- matrix(y,ncol=1)
	
	yvar <- mf$"(variance)"
		
 	res <- metareg(y,x,yvar,threshold=threshold, maxiter=maxiter, npermut=npermut)
	
	res$call <- call
	
	class(res) <- "metatest"
	
	return(res)
}

