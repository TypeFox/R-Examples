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

