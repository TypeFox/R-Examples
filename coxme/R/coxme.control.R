#
# Gather all of the control parameters for coxme into one spot
#
coxme.control <- function(eps=1e-8, 
                          toler.chol = .Machine$double.eps ^ .75, 
                          iter.max =20,
			  inner.iter=Quote(max(4, fit0$iter+1)),
			  sparse.calc=NULL,
                          optpar=list(method='BFGS', 
                                      control=list(reltol=1e-5)),
                          refine.df = 4, refine.detail=FALSE,
                          refine.method="control",
                          sparse=c(50, .02), 
                          varinit = c(.02, .1, .4, .8)^2,
                          corinit = c(0, .3)) {
    if (iter.max <0) stop("Invalid value for iterations")
    if (inner.iter<1) stop("Invalid value for inner iterations")
    if (eps <=0) stop ("Invalid convergence criteria")
    if (eps <= toler.chol) 
	    warning("For numerical accuracy, tolerance should be < eps")
    if (optpar$control$reltol <= eps)
        warning(paste("For numerical accuracy, eps (tolerance for an inner",
                      "loop) should be < the relative tol of optim"))
    if (!is.null(sparse.calc)) {
        if (sparse.calc !=0 && sparse.calc !=1)
            stop("Invalid value for sparse.calc option")
        }
    if (length(varinit)<1 || !is.numeric(varinit))
        stop("varinit must be a vector of numeric values")
    if (length(corinit)<1 || !is.numeric(corinit))
        stop("corinit must be a vector of numeric values")
    if (any(varinit <=0)) stop ("varinit values must be >0")
    if (any(corinit <0))  stop("corinit values must be >=0")

    list(eps=eps, toler.chol=toler.chol, iter.max=iter.max,
	 inner.iter=inner.iter, sparse.calc=sparse.calc,
         optpar=optpar, refine.df=refine.df, refine.detail=refine.detail,
         refine.method=refine.method,
         sparse=sparse, varinit=varinit, corinit=corinit)
    }
