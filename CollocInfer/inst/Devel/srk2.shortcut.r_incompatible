Smooth.SRK2 <- function(fn, data, times, pars, coefs=NULL, basisvals=NULL,
                       lambda,fd.obj=NULL,more=NULL,weights=NULL,
	                     quadrature=NULL, in.meth='nlminb', control.in=list(), eps=1e-6,
                       pos=0, discrete=0, names=NULL)
{
      
    dims = dim(data)

    profile.obj = SRK2.setup(pars, coefs, fn, basisvals, lambda,fd.obj, more,
                            data, weights, times, quadrature, eps=1e-6, pos,
                            discrete, names=names)

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


Profile.SRK2 <- function(fn,data,times,pars,coefs=NULL,basisvals=NULL,lambda,
                        fd.obj=NULL,more=NULL,weights=NULL,quadrature=NULL,
                        in.meth='nlminb',out.meth='nls',
                        control.in=list(),control.out=list(),eps=1e-6,
                        active=NULL,pos=0,discrete=0,names=NULL)
{
    if(is.null(active)){ active = 1:length(pars) }

    profile.obj = SRK2.setup(pars,coefs,fn,basisvals,lambda,fd.obj,more,
      data,weights,times,quadrature,eps=1e-6,pos,discrete,names)

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
    
    if(out.meth == "ProfileGN"){
      res=Profile.GausNewt(pars=pars,times=times,data=data,coefs=ncoefs,
		    lik=lik,proc=proc,in.meth=in.meth,control.in=control.in,
		    active=active,control=control.out)
      apars = res$pars[active]

      ncoefs = res$in.res$coefs
      g = res$in.res$df
      resid = res$in.res$f
    }
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

SRK2.setup = function(pars, coefs=NULL, fn, basisvals=NULL, lambda, fd.obj=NULL,
                     more=NULL, data=NULL, weights=NULL, times=NULL,
                     quadrature=NULL, eps=1e-6, pos=0, discrete=0, names=NULL)
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
    
    if(pos==0) 
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
    else{
      stop('fn must be either a list of functions or a function')
    }
    
    if(pos==0) proc$more = procmore
    else { proc$more = make.logtrans()
           proc$more$more = procmore}
           
    
    proc$more$names = colnames
    proc$more$parnames = names(pars)


    if(is.basis(basisvals)){
      if(is.null(times)){
        stop(paste('if basisvals is is a basis object,',
                   ' you must specify the observation times'))
      }

      lik$bvals = eval.basis(times,basisvals)
           
      if(is.null(quadrature) | is.null(quadrature$qpts)){
        knots = c(basisvals$rangeval[1],basisvals$params,basisvals$rangeval[2])
        qpts = knots
      }
      else{
        qpts = quadrature$qpts
      }

        proc$bvals = list(bvals = Matrix(eval.basis(qpts,basisvals,0),sparse=TRUE),
                          I = SRK2fns()$SRK2indeces(length(qpts)))
        proc$more$weights = matrix(1,length(qpts),ncol(coefs))
        proc$more$qpts = qpts
    }
    else{                                   # quadrature is ignored if basisvals
                                            # is not a basis object
      if(is.matrix(basisvals) | is.null(basisvals)){
        if(is.null(basisvals)){ basisvals = Diagonal(nrow(coefs)) }
        lik$bvals = basisvals
        proc$bvals = list(bvals  = bvals,
                          I = SRK2fns()$SRK2indeces(nrow(coefs)))
        proc$more$weights = matrix(1,nrow(basisvals)-1,ncol(coefs))
        proc$more$qpts = times[1:(length(times)-1)]
      }                                    
      else{                                      
        lik$bvals = basisvals$bvals.obs
  
        proc$bvals =  list(bvals=basisvals$bvals,
                          I = SRK2fns()$SRK2indeces(nrow(basisvals$bvals)))
        proc$more$qpts = basisvals$qpts
        
        if(!is.null(basisvals$qwts))
          proc$more$weights = matrix(basisvals$qwts, 
                                     length(basisvals$qwts),ncol(coefs))
        else
          proc$more$weights = matrix(1/,length(proc$more$qpts),ncol(coefs))    
      }
    }
 
    if( is.null(weights) ){
      lik$more$weights = matrix(1,ncol(coefs))
    }
    else{
      lik$more$weights = matrix(weights,nrow(lik$bvals),ncol(coefs))
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

    if(length(lambda)==1){ lambda = rep(lambda,ncol(coefs)) }

#    print(length(proc$more$weights))
    proc$more$weights = proc$more$weights%*%diag(lambda)

    return( list(lik=lik,proc=proc,coefs=coefs,data=data,times=times) )

}
