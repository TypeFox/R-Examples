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

###############################################################################


Profile.LS <- function(fn,data,times,pars,coefs=NULL,basisvals=NULL,lambda,
                        fd.obj=NULL,more=NULL,weights=NULL,quadrature=NULL,
                        in.meth='nlminb',out.meth='nls',
                        control.in=list(),control.out=list(),eps=1e-6,
                        active=NULL,posproc=FALSE,poslik=FALSE,discrete=FALSE,names=NULL,sparse=FALSE)
{
#    browser()
    if(is.null(active)){ active = 1:length(pars) }

    profile.obj = LS.setup(pars,coefs,fn,basisvals,lambda,fd.obj,more,
      data,weights,times,quadrature,eps=1e-6,posproc,poslik,discrete,names,sparse)

    dims = dim(data)

    lik   = profile.obj$lik
    proc  = profile.obj$proc
    coefs = profile.obj$coefs
    data  = profile.obj$data
    times = profile.obj$times
   
    if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
    if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
    if(file.exists('counter.tmp')){file.remove('counter.tmp')}
    
    Ires = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = Ires$coefs

    write.table(ncoefs,file='optcoefs.tmp',col.names=FALSE,row.names=FALSE)
    write.table(ncoefs,file='curcoefs.tmp',col.names=FALSE,row.names=FALSE)
    
    apars = pars[active]
    aparamnames = names(apars)
    
    ################  Gauss-Newton optimization  #########################
    if(out.meth == "ProfileGN"){
      res=Profile.GausNewt(pars=pars, times=times, data=data, coefs=ncoefs,
		               lik=lik, proc=proc, in.meth=in.meth,
                           control.in=control.in,
		               active=active, control=control.out)
      apars = res$pars[active]

      ncoefs = res$in.res$coefs
      g = res$in.res$df
      resid = res$in.res$f
    }

    ################  nls optimization  #########################

    if(out.meth == "nls"){
      if(is.null(control.out$trace)){control.out$trace=TRUE}
      if(is.null(control.out$maxiter)){control.out$maxiter=100}
      if(is.null(control.out$tol)){control.out$tol=1e-8}
      if(is.null(control.out$printEval)){control.out$printEval=TRUE}
      if(is.null(control.out$warnOnly)){control.out$warnOnly=TRUE}
      res = nls(~ProfileSSE(pars, allpars, times, data, coefs, lik, proc,
                            in.meth, control.in, active),
        data = list(allpars=pars, times=times, data=data, coefs=ncoefs,
                    lik=lik, proc=proc,
        in.meth=in.meth,control.in=control.in,active=active),
        start = list(pars=pars[active]),trace=control.out$trace,control=control.out)   
      apars = res$m$getPars()

      g = res$m$gradient()
      resid = res$m$resid()
      if(file.exists('curcoefs.tmp'))
      	 ncoefs = as.matrix(read.table(file='curcoefs.tmp'))
      else 
         ncoefs = coefs
    }

    names(apars) = aparamnames

    pars[active] = apars

    ncoefs = as.matrix(ncoefs)
    if(!is.null(proc$more$names)){ colnames(ncoefs) = proc$more$names }

   
     if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
     if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
     if(file.exists('counter.tmp')){file.remove('counter.tmp')}
    
    if(!is.null(fd.obj)){
      ncoefs = array(ncoefs,c(nrow(ncoefs)/dims[2],dims[2],dims[3]))
      fd.obj = fd(ncoefs,fd.obj$basis)
      return( list(pars=pars,fd=fd.obj,lik=lik,proc=proc,outer.result=res) )
    }
    else{
      return( list(pars=pars, coefs=ncoefs, lik=lik, proc=proc, outer.result=res,
                   data=data, times=times) )
    }    
}

###############################################################################

LS.setup = function(pars, coefs=NULL, fn, basisvals=NULL, lambda, fd.obj=NULL,
                     more=NULL, data=NULL, weights=NULL, times=NULL,
                     quadrature=NULL, eps=1e-6, posproc=FALSE, poslik = FALSE, 
                     discrete=FALSE, names=NULL, sparse=FALSE)
{
    colnames = names

    if(!is.null(fd.obj)){            # If an fd object is provided, it overrides
      basisvals = fd.obj$basis       # the basis and function values
      
      if(!is.null(fd.obj$coefs)){
        coefs = fd.obj$coefs
      }
      if(!is.null(fd.obj$fdnames) & is.null(colnames)){
        colnames = fd.obj$fdnames[[length(fd.obj$fdnames)]]
      }
    }
        
    lik = make.SSElik()
    
    if(!poslik) 
       lik$more = make.id()
    else
      lik$more = make.exp()
    
    
    if(length(dim(coefs))>2){
      if(is.null(colnames)){
        colnames = dimnames(coefs)[[3]]
      }
      nrep = dim(coefs)[2]
      coefs = matrix(coefs,dim(coefs)[1]*dim(coefs)[2],dim(coefs)[3])
    }
    else{
      nrep = 1
      if(is.null(colnames)){
        colnames = colnames(coefs)
      }
    }


    proc = make.SSEproc()
    
    if(is.list(fn)){
      procmore = fn
      procmore$more = more
    }
    else if(is.function(fn)){
      procmore = make.findif.ode()
      procmore$more$fn = fn
      procmore$more$more = more
      procmore$more$eps = eps
    }          
    else if(inherits(fn,'pomp')){
      procmore = make.findif.ode()
      procmore$more$fn = pomp.skeleton
      procmore$more$eps = eps
      procmore$more$more =  list(pomp.obj = fn)
    }
    else{
      stop('fn must be either a list of functions or a function')
    }
    
    if(!posproc) proc$more = procmore
    else { proc$more = make.logtrans()
           proc$more$more = procmore}
           
    
    proc$more$names = colnames
    proc$more$parnames = names(pars)


    if(is.basis(basisvals)){
      if(is.null(times)){
        stop(paste('if basisvals is is a basis object,',
                   ' you must specify the observation times'))
      }

      if(sparse){
        lik$bvals = spam(diag(rep(1,nrep)) %x% 
                   eval.basis(times,basisvals))
      } else{
        lik$bvals = diag(rep(1,nrep)) %x% eval.basis(times,basisvals)
      }      
           
      if(is.null(quadrature) | is.null(quadrature$qpts)){
        knots = c(basisvals$rangeval[1],basisvals$params,basisvals$rangeval[2])
        qpts = knots[-length(knots)] + diff(knots)/2
      }
      else{
        qpts = quadrature$qpts
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
        proc$more$weights = matrix(1/length(qpts),length(qpts)*nrep,ncol(coefs))
        #proc$more$weights = 1/length(qpts)
        proc$more$qpts = qpts
      }
      else{
       len = length(times)
       basis = eval.basis(times,basisvals,0)
       if(sparse){
         proc$bvals = list(bvals = spam(basis[1:(len-1),]),
                           dbvals= spam(basis[2:len,]))
       } else{
         proc$bvals = list(bvals = basis[1:(len-1),],
                           dbvals= basis[2:len,])       
       }
       proc$more$weights = matrix(1/(len-1),(len-1)*nrep,ncol(coefs))
       #proc$more$weights = 1/(len-1)
       proc$more$qpts = times[1:(len-1)]
       }
      
      
    }
    else{                                   # quadrature is ignored if basisvals
                                            # is not a basis object
      if(discrete & (is.matrix(basisvals) | is.null(basisvals))){
        if(is.null(basisvals)){ basisvals = Diagonal(nrow(coefs)) }
        if(sparse){
          lik$bvals = spam(diag(rep(1,nrep)) %x%basisvals) 
          proc$bvals = list(bvals  = spam(diag(rep(1,nrep)) %x% basisvals[1:(nrow(basisvals)-1),]),
                            dbvals = spam(diag(rep(1,nrep)) %x% basisvals[2:nrow(basisvals),]))        
        }else{
          lik$bvals = diag(rep(1,nrep)) %x%basisvals  
          proc$bvals = list(bvals  = diag(rep(1,nrep)) %x% basisvals[1:(nrow(basisvals)-1),],
                            dbvals =diag(rep(1,nrep)) %x% basisvals[2:nrow(basisvals),])
        }
        proc$more$weights = matrix(1/(nrow(basisvals)-1),
                                   nrow(basisvals)-1,ncol(coefs))
        #proc$more$weights = 1/(nrow(basisvals)-1)
        proc$more$qpts = times[1:(length(times)-1)]
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
        
        if(!is.null(basisvals$qwts))
          proc$more$weights = matrix(basisvals$qwts, 
                                     length(basisvals$qwts)*nrep,ncol(coefs))
        else
          #proc$more$weights = 1/length(proc$more$qpts)
          proc$more$weights = matrix(1/length(proc$more$qpts),
                                       length(proc$more$qpts)*nrep,ncol(coefs))    
      }
    }
 
    if( is.null(weights) ){
      lik$more$weights = matrix(1,ncol(coefs))
    }
    else{
      lik$more$weights = matrix(weights,nrow(lik$bvals)*nrep,ncol(coefs))
    }

    if(!is.null(data)){
      if(length(dim(data))==2){
        if(nrep>1){stop('data dimensions must match coefficient dimensions')}
        if(dim(data)[1] != length(times) | dim(data)[2]!= dim(coefs)[2]){
        stop('data dimensions, times and coefficient dimensions do not match')}
      }
      if(length(dim(data))==3){
         if(dim(data)[2] != nrep | dim(data)[3]!=dim(coefs)[2] | 
            dim(data)[1]!=length(times)){
        stop('data dimensions, times and coefficient dimensions do not match')}
        data = matrix(data,dim(data)[1]*dim(data)[2],dim(data)[3])
        times = rep(times,nrep)
     }
    }

 #   if(length(lambda)==1){ lambda = as.matrix(rep(lambda,ncol(coefs))) }

    if(length(lambda) > 1){ proc$more$weights = proc$more$weights%*%diag(lambda) }
    else{ proc$more$weights = as.numeric(lambda)*proc$more$weights }

    return( list(lik=lik,proc=proc,coefs=coefs,data=data,times=times) )

}
