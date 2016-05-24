#  This file defines three functions:

#  Profile.multinorm:  Optimizing parameter values using outer optimization
#  Smooth.multinorm:   Smoothing data using the inner coefficient optimization
#  multinorm.setup:    Set up multivariate normal negative log density functions

################################################################################

Profile.multinorm <- function(fn,data,times,pars,coefs=NULL,basisvals=NULL,var=c(1,0.01),
                              fd.obj=NULL,more=NULL,quadrature=NULL,
                              in.meth='nlminb',out.meth='optim',
                              control.in,control.out,eps=1e-6,active=NULL,pos=0,
                              discrete=0)
{
    if(is.null(active)){ active = 1:length(pars) }

    profile.obj = multinorm.setup(pars,coefs,fn,basisvals,var,fd.obj,more,
                                  times,quadrature,eps=1e-6,pos=pos,
                                  discrete=discrete)
    
    lik   = profile.obj$lik
    proc  = profile.obj$proc
    coefs = profile.obj$coefs
         
    Ires = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    
    apars = pars[active]
    aparamnames = names(apars)
    
    res = outeropt(data,times,pars,Ires$coefs,lik,proc,in.meth,out.meth,
                   control.in,control.out,active)
    
    apars = res$pars[active]
    names(apars) = aparamnames
    ncoefs = as.matrix(res$coefs)
    
    pars[active] = apars
    
    if(!is.null(fd.obj)){
      fd.obj = fd(ncoefs,fd.obj$basis)
      return( list(pars=pars,fd=fd.obj,lik=lik,proc=proc,res=res$res) )
    }
    else{
      return( list(pars=pars,coefs=ncoefs,lik=lik,proc=proc,res$res) )
    }
}

################################################################################

Smooth.multinorm <- function(fn,data,times,pars,coefs=NULL,basisvals=NULL,var=c(1,0.01),
                             fd.obj=NULL,more=NULL,quadrature=NULL,
                             in.meth='nlminb',control.in,eps=1e-6,pos=0,discrete=0)
{

    profile.obj = multinorm.setup(pars,coefs,fn,basisvals,var,fd.obj,more,
                                  times,quadrature,eps=1e-6,pos=pos,discrete=discrete)

    lik  = profile.obj$lik
    proc = profile.obj$proc
    coefs = profile.obj$coefs
    
    dims = dim(data)

    if(length(dims)>2){
      data = matrix(data,dim(data)[1]*dim(data)[2],dim(data)[3])
      dims[1] = length(coefs)/(dims[2]*dims[3])
    }
    else{
      dims = dim(coefs)
    }

    Ires = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = Ires$coefs
    Ires = Ires$res
    ncoefs = array(ncoefs,dims)

    if(!is.null(fd.obj)){
      fd.obj = fd(ncoefs,fd.obj$basis)
      return( list(fd=fd.obj,lik=lik,proc=proc,res=Ires) )
    }
    else{
      return( list(coefs=ncoefs,lik=lik,proc=proc,res=Ires) )
    }


}

################################################################################

multinorm.setup = function(pars,coefs=NULL,fn,basisvals=NULL,var=c(1,0.01),fd.obj=NULL,
                           more=NULL,times=NULL,quadrature=NULL,eps=1e-6,pos=0,discrete=0)
{

    if(!is.null(fd.obj)){                 # If an fd object is provided, it overrides
      basisvals = fd.obj$basis            # the basis and function values

      if(!is.null(fd.obj$coefs)){
        coefs = fd.obj$coefs
      }

      colnames = fd.obj$fdnames[[3]]
    }
    else{
      colnames = NULL
    }


    lik = make.multinorm()
    
    if(pos==0){ lik$more = c(make.id(),make.cvar())}
    else { lik$more = c(make.exp(),make.cvar())}
    
    lik$more$f.more = NULL
    lik$more$v.more= list(mat=var[1]*diag(rep(1,2)),sub=matrix(0,0,3))
    
    if(length(dim(coefs))>2){
      if(is.null(colnames)){
        colnames = dimnames(coefs)[[3]]
      }
      nrep = dim(coefs)[2]
      coefs = matrix(coefs,dim(coefs)[1]*dim(coefs)[2],dim(coefs)[3])
    }
    else{
      nrep = 1
      colnames = colnames(coefs)
    }

    if(pos==0){
    	if(discrete==0) proc = make.Cproc()
    	else 
    	 proc = make.Dproc()
    	}
    else{
    	if(discrete==0) proc = make.exp.Cproc()
    	else 
    	 proc = make.exp.Dproc()
    }
    
    proc$more = make.multinorm()

    if(is.list(fn)){
      proc$more$more = c(fn,make.cvar())
      proc$more$more$f.more = NULL
      proc$more$more$v.more = list(mat=var[2]*diag(rep(1,2)),sub=matrix(0,0,3))
      proc$more$more$more = more
    }
    else if(is.function(fn)){
      proc$more$more = c(make.findif.ode(),make.cvar())
      proc$more$more$f.more$eps = eps
      proc$more$more$f.more$fn = fn
      proc$more$more$more = more
    }
    else{
      stop('fn must be either a list of functions or a function')
    }

    proc$more$names = colnames
    proc$more$parnames = names(pars)


    if(is.basis(basisvals)){
      if(is.null(times)){
        stop('if basisvals is is a basis object, you must specify the observation times')
      }

       lik$bvals = Matrix(diag(rep(1,nrep)) %x% eval.basis(times,basisvals),sparse=TRUE)
       
      if(is.null(quadrature) | is.null(quadrature$qpts)){
        knots = c(basisvals$rangeval[1],basisvals$params,basisvals$rangeval[2])
        qpts = knots[-length(knots)] + diff(knots)
      
      }
      else{
        qpts = quadrature$qpts
      }

      proc$bvals = list()
     
      if(discrete==0){
        proc$bvals$bvals  = Matrix(diag(rep(1,nrep))%x%eval.basis(qpts,basisvals,0),sparse=TRUE)
        proc$bvals$dbvals = Matrix(diag(rep(1,nrep))%x%eval.basis(qpts,basisvals,1),sparse=TRUE)
        proc$more$qpts = qpts
      }
      else{
       len = length(times)
       proc$bvals = Matrix(diag(rep(1,nrep))%x%eval.basis(times,basisvals,0),sparse=TRUE)
       proc$more$qpts = times[1:(len-1)]
      }
    }
    else{    
      if(discrete & (is.matrix(basisvals) | is.null(basisvals))){
        if(is.null(basisvals)){ basisvals = Diagonal(nrow(coefs)) }
        lik$bvals = basisvals
        proc$bvals = basisvals
        proc$more$qpts = times[1:(length(times)-1)]
      }                                    
      else{                                  
        lik$bvals = basisvals$bvals.obs
  
        proc$bvals =  list(bvals=basisvals$bvals,dbvals=basisvals$dbvals)
        proc$more$qpts = basisvals$qpts
      } 
       
    }

    return( list(lik=lik,proc=proc,coefs=coefs) )

}
