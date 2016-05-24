Smooth.LS <- function(fn, data, times, pars, coefs=NULL, basisvals=NULL,
                       lambda,fd.obj=NULL,more=NULL,weights=NULL,
	                     quadrature=NULL, in.meth='nlminb', control.in=list(), eps=1e-6,
                       posproc=FALSE, poslik=FALSE, discrete=FALSE, names=NULL, sparse=FALSE)
{
      
    dims = dim(data)

    profile.obj = LS.setup(pars, coefs, fn, basisvals, lambda,fd.obj, more,
                            data, weights, times, quadrature, eps=1e-6, posproc,
                            poslik, discrete, names=names, sparse=sparse)

    lik   = profile.obj$lik
    proc  = profile.obj$proc
    coefs = profile.obj$coefs
    data  = profile.obj$data
    times = profile.obj$times

    Ires   = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = Ires$coefs
    Ires   = Ires$res
    ncoefs = as.matrix(ncoefs)
    if(!is.null(proc$more$names)){ colnames(ncoefs) = proc$more$names }

    if(!is.null(fd.obj)){
      if(length(dims)>2){
        ncoefs = array(ncoefs,c(length(ncoefs)/(dims[2]*dims[3]),dims[2],dims[3]))
      } else{
         ncoefs = array(ncoefs,c(length(ncoefs)/dims[2],dims[2]))
      }
      fd.obj = fd(ncoefs,fd.obj$basis)
      return( list(fd=fd.obj,lik=lik,proc=proc,inner.result=Ires) )
    }
    else{
      return( list(coefs=ncoefs,lik=lik,proc=proc,inner.result=Ires,data=data,
                   times=times) )
    }
}


