ldr <-
function(X, y=NULL, fy=NULL, Sigmas=NULL, ns=NULL, numdir=NULL, nslices=NULL, 
model=c("core", "lad", "pfc"), numdir.test=FALSE, ...)
{
	if (model=="pfc")
	{	
		if (is.null(fy)){stop("fy is not provided"); return()}
 
 		return(invisible(pfc(X=X, y=y, fy=fy, numdir=numdir, numdir.test=numdir.test, ...)))
	}

	if (model=="lad") 
	{
		if (is.null(y)){stop("The response is needed"); return()}
	
		return(invisible(lad(X=X, y=y, numdir=numdir, nslices=nslices, numdir.test=numdir.test,...)))
	}

	if (model=="core") 
	{
		if (!is.null(Sigmas) & !is.null(ns)) return(invisible(core(Sigmas=Sigmas, ns=ns, numdir=numdir, numdir.test=numdir.test,...)))
        
        	return(invisible(core(X=X, y=y, numdir=numdir, numdir.test=numdir.test,...)))
	}
}
