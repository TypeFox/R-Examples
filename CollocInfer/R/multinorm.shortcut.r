#  This file defines three functions:

#  Profile.multinorm:  Optimizing parameter values using outer optimization
#  Smooth.multinorm:   Smoothing data using the inner coefficient optimization
#  multinorm.setup:    Set up multivariate normal negative log density functions

################################################################################

Profile.multinorm <- function(fn,data,times,pars,coefs=NULL,basisvals=NULL,var=c(1,0.01),
                        fd.obj=NULL,more=NULL,quadrature=NULL,
                        in.meth='nlminb',out.meth='optim',
                        control.in=list(),control.out=list(),eps=1e-6,
                        active=NULL,posproc=FALSE,poslik=FALSE,discrete=FALSE,names=NULL,sparse=FALSE)
{

    dims = dim(data)
    if(is.null(active)){ active = 1:length(pars) }


    profile.obj = multinorm.setup(pars,coefs,fn,basisvals,var,fd.obj,more,
          data,times,quadrature,eps=1e-6,posproc=posproc,poslik=poslik,
          discrete=discrete,names=names,sparse=sparse)
    
    lik  = profile.obj$lik
    proc = profile.obj$proc
    coefs = profile.obj$coefs
    data  = profile.obj$data
    times = profile.obj$times
    
    Ires = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    
    apars = pars[active]
    aparamnames = names(apars)
    
    res = outeropt(data,times,pars,Ires$coefs,lik,proc,in.meth,out.meth,control.in,control.out,active)
    
    apars = res$pars[active]
    names(apars) = aparamnames
    ncoefs = as.matrix(res$coefs)
    
    pars[active] = apars

    if(!is.null(fd.obj)){
      if(length(dims)>2){
        ncoefs = array(ncoefs,c(length(ncoefs)/(dims[2]*dims[3]),dims[2],dims[3]))
      } else{
         ncoefs = array(ncoefs,c(length(ncoefs)/dims[2],dims[2]))
      }
      fd.obj = fd(ncoefs,fd.obj$basis)
      return( list(pars=pars,fd=fd.obj,lik=lik,proc=proc,outer.result=res$outer.result,data=data,times=times) )
    }
    else{
      return( list(pars=pars,coefs=ncoefs,lik=lik,proc=proc,outer.result=res$outer.result,data=data,times=times) )
    }
}

################################################################################

Smooth.multinorm <- function(fn,data,times,pars,coefs=NULL,basisvals=NULL,var=c(1,0.01),
                        fd.obj=NULL,more=NULL,quadrature=NULL,
                        in.meth='nlminb',control.in=list(),
                        eps=1e-6,posproc=FALSE,poslik=FALSE,discrete=FALSE,names=NULL,sparse=FALSE)
{

    dims = dim(data)
    profile.obj = multinorm.setup(pars,coefs,fn,basisvals,var,fd.obj,more,
        data,times,quadrature,eps=1e-6,posproc=posproc,poslik=poslik,
        discrete=discrete,names=names,sparse=sparse)

    lik  = profile.obj$lik
    proc = profile.obj$proc
    coefs = profile.obj$coefs
    data  = profile.obj$data
    times = profile.obj$times
    
    dims = dim(data)

    Ires = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = Ires$coefs
    Ires = Ires$res
    ncoefs = array(ncoefs,dims)

    if(!is.null(fd.obj)){
      ncoefs = array(ncoefs,c(nrow(ncoefs)/dims[2],dims[2],dims[3]))
      fd.obj = fd(ncoefs,fd.obj$basis)
      return( list(fd=fd.obj,lik=lik,proc=proc,res=Ires) )
    }
    else{
      return( list(coefs=ncoefs,lik=lik,proc=proc,inner.result=Ires,data=data,times=times) )
    }


}

################################################################################

multinorm.setup = function(pars,coefs=NULL,fn,basisvals=NULL,var=c(1,0.01),fd.obj=NULL,
  more=NULL,data=NULL,times=NULL,quadrature=NULL,eps=1e-6,posproc=FALSE,poslik=FALSE,
  discrete=FALSE,names=NULL,sparse=FALSE)
{

    if(!is.null(data) & length(dim(data))>2){
        data = matrix(data,dim(data)[1]*dim(data)[2],dim(data)[3])
    }
    
    colnames = names
    if(!is.null(fd.obj)){                 # If an fd object is provided, it overrides
      basisvals = fd.obj$basis            # the basis and function values

      if(!is.null(fd.obj$coefs)){
        coefs = fd.obj$coefs
      }
      if(!is.null(fd.obj$fdnames) & is.null(colnames)){
        colnames = fd.obj$fdnames[[length(fd.obj$fdnames)]]
      }
    }


    lik = make.multinorm()
    
    if(!poslik){ lik$more = c(make.id(),make.cvar())}
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

    if(!posproc){
    	if(!discrete) proc = make.Cproc()
    	else 
    	 proc = make.Dproc()
    	}
    else{
    	if(!discrete) proc = make.exp.Cproc()
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
    else if(inherits(fn,'pomp')){
      proc$more$more = c(make.findif.ode(),make.cvar())
      proc$more$more$f.more$fn = pomp.skeleton
      proc$more$more$f.more$eps = eps
      proc$more$more$f.more$more =  list(pomp.obj = fn)
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

      if(sparse){
        lik$bvals = spam(diag(rep(1,nrep)) %x% 
                   eval.basis(times,basisvals))
      } else{
        lik$bvals = diag(rep(1,nrep)) %x% eval.basis(times,basisvals)
      }      
             
      if(is.null(quadrature) | is.null(quadrature$qpts)){
        knots = c(basisvals$rangeval[1],basisvals$params,basisvals$rangeval[2])
        qpts = c(knots[1],knots[-length(knots)] + diff(knots)/2,knots[length(knots)])
        weights = rep(1,length(qpts))
      
      }
      else{
        qpts = quadrature$qpts
        weights = quadrature$weights
        if(is.null(weights)){ weights = rep(1,length(qpts)) }
      }

      proc$bvals = list()
     
      if(!discrete){
        if(sparse){
          proc$bvals$bvals  = spam(diag(rep(1,nrep)) %x% 
                                     eval.basis(qpts,basisvals,0))
          proc$bvals$dbvals = spam(diag(rep(1,nrep)) %x%
                                     eval.basis(qpts,basisvals,1))
        }else{
          proc$bvals$bvals  = diag(rep(1,nrep)) %x% eval.basis(qpts,basisvals,0)
          proc$bvals$dbvals = diag(rep(1,nrep)) %x% eval.basis(qpts,basisvals,1)        
        }
        proc$more$qpts = rep(qpts,nrep)
        proc$mroe$weights = rep(weights,nrep)
      }
      else{
       len = length(times)
       if(sparse){
         proc$bvals = list(bvals = spam(basisvals[1:(len-1),]),
                           dbvals= spam(basisvals[2:len,]))
       } else{
         proc$bvals = list(bvals = basisvals[1:(len-1),],
                           dbvals= basisvals[2:len,])
       }
       proc$more$qpts = rep(times[1:(len-1)],nrep)
      }
    }
    else{    
      if(discrete & (is.matrix(basisvals) | is.null(basisvals))){
        if(is.null(basisvals)){ basisvals = Diagonal(nrow(coefs)) }
        if(sparse){
          lik$bvals = spam(diag(rep(1,nrep))%x%basisvals)
          proc$bvals = spam(diag(rep(1,nrep))%x%basisvals)
        }else{
          lik$bvals = diag(rep(1,nrep))%x%basisvals
          proc$bvals = diag(rep(1,nrep))%x%basisvals      
        }
        proc$more$qpts = rep(times[1:(length(times)-1)],nrep)
      }                                    
      else{                                  
        if(sparse){                                 
          lik$bvals = spam(diag(rep(1,nrep))%x%basisvals$bvals.obs)
    
          proc$bvals =  list(bvals=spam(diag(rep(1,nrep)) %x% 
                                          basisvals$bvals),
                            dbvals=spam(diag(rep(1,nrep)) %x%
                                          basisvals$dbvals))
        } else{
          lik$bvals = diag(rep(1,nrep))%x%basisvals$bvals.obs
    
          proc$bvals =  list(bvals=diag(rep(1,nrep)) %x% basisvals$bvals,
                            dbvals= diag(rep(1,nrep)) %x%basisvals$dbvals)        
        }
        proc$more$qpts = rep(basisvals$qpts,nrep)
        proc$more$weights = rep(basisvals$weights,nrep)
      } 
       
    }
    if(is.null(proc$more$weights)){ proc$more$weights = rep(1,length(proc$more$qpts)) }
    
    if(!is.null(data)){
      if(length(dim(data))==2){
        if(nrep>1){stop('data dimensions must match coefficient dimensions')}
        if(dim(data)[1] != length(times) | dim(data)[2]!= dim(coefs)[2]){stop('data dimensions, times and coefficient dimensions do not match')}
      }
      if(length(dim(data))==3){
         if(dim(data)[2] != nrep | dim(data)[3]!=dim(coefs)[2] | dim(data)[1]!=length(times)){
                      stop('data dimensions, times and coefficient dimensions do not match')}
        data = matrix(data,dim(data)[1]*dim(data)[2],dim(data)[3])
        times = rep(times,nrep)
     }
    }
    
    

    return( list(lik=lik,proc=proc,coefs=coefs,data=data,times=times) )

}
