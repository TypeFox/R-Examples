#  This file defines four functions:
#  inneropt:  optimization with respect to coefficients given parameter values
#  outeropt:  optimization of parameters
#  FitMatchOpt:  
#  ParsMatchOpt:

##############################################################################

inneropt <- function(data, times, pars, coefs, lik, proc, in.meth='nlminb', 
                     control.in=NULL)
{

#  The multivariate data in argment data are observations of one or more
#  functional variables or processes at time points in argument time.

#  The data are smoothed using differential operators.  These operators are 
#  defined by the functions in argument fn.  

#  The function optimizes the coefficients of the basis expansions for the
#  representations of the variables or processes in what is called the
#  inner optimization loop.

#  The fitting criterion for both the data and the differential equation
#  is user-defined in arguments lik and proc, respectively.

#  The optimization method for the inner loop is defined by argument in.meth.

#  The function returns the optimized coefficients and the residual values.

  if (in.meth=="SplineEst"){
    #  optimization by Newton-Raphson using function SplineEst.NewtRaph
        res = SplineEst.NewtRaph(coefs,times,data,lik,proc,pars,control.in)
        ncoefs = matrix(res$coefs,ncol(lik$bvals),length(res$coefs)/ncol(lik$bvals))
    }
    else if(in.meth=="optim"){
    #  optimization using function optim
        res = optim(coefs,SplineCoefsErr,gr=SplineCoefsDC,hessian=T,
            control=control.in,times=times,data=data,lik=lik,proc=proc,pars=pars,method="BFGS")

        ncoefs = matrix(res$par,ncol(lik$bvals),length(res$par)/ncol(lik$bvals))
    }    
    else if(in.meth=="nlminb"){
    #  optimization using function nlminb
        res = nlminb(coefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=SplineCoefsDC2,
            control=control.in,times=times,data=data,lik=lik,proc=proc,pars=pars)

        ncoefs = matrix(res$par,ncol(lik$bvals),length(res$par)/ncol(lik$bvals))
    }
    else if(in.meth=="maxNR"){
    #  optimizaton using function maxNR
        if(is.null(control.in$print.level)){control.in$print.level = 1}
        if(is.null(control.in$iterlim)){control.in$iterlim = 1000}    
        res = maxNR(SplineCoefsErr,coefs,times=times,data=data,lik=lik,proc=proc,pars=pars,sgn=-1,
            grad=SplineCoefsDC,hess=SplineCoefsDC2,print.level=control.in$print.level,
            iterlim = control.in$iterlim)

        ncoefs = matrix(res$estimate,ncol(lik$bvals),length(res$estimate)/ncol(lik$bvals))
    }
    else if(in.meth=='trust'){
    #  optimization using functon trust
        if(is.null(control.in$rinit)){ control.in$rinit = 1}
        if(is.null(control.in$rmax)){ control.in$rinit = 100}
        if(is.null(control.in$parscale)){ control.in$parscale = abs(coefs)}
        if(is.null(control.in$iterlim)){ control.in$iterlim = 100}
        
        res = trust(SplineCoefsList,coefs,rinit=control.in$rinit,rmax=control.in$rmax,parscale=control.in$parscale,iterlim=control.in$iterlim,
          times=time,data=data,lik=lik,proc=proc,pars=pars,sgn=1)
    
        ncoefs = matrix(res$argument,dim(coefs))
    }

    if(!is.null(proc$more$names)){ colnames(ncoefs) = proc$more$names }
  
    return(list(coefs=ncoefs,res=res))
}

##############################################################################

outeropt <- function(data, times, pars, coefs, lik, proc,
                     in.meth='nlminb', out.meth='nlminb',
                     control.in=NULL, control.out=NULL,
                     active=1:length(pars))
{

#  The multivariate data in argment data are observations of one or more
#  functional variables or processes at time points in argument times.

#  The data are smoothed using differential operators.  

#  The operators typically depend on a number of parameters whose values
#  are in argument pars.  The function optimizes these parameter values.

#  The fitting criterion for both the data and the differential equation
#  is user-defined in arguments lik and proc, respectively.

#  The function optimizes the coefficients of the basis expansions for the
#  representations of the variables or processes in what is called the
#  inner optimization loop, implemented by function Inneropt.  It optimizes the
#  parameter values in the outer optimization loop.

#  The optimization method for the inner loop is defined by argument in.meth.
#  The optimization method for the outer loop is defined by argument out.meth.

#  The function returns the optimized parameter values, a functional data object defining 
#  the smoothing functions, the likelihood and process functions, and the
#  residual values.

    if(file.exists('curcoefs.tmp')) file.remove('curcoefs.tmp')
    if(file.exists('optcoefs.tmp')) file.remove('optcoefs.tmp')
    if(file.exists('counter.tmp'))  file.remove('counter.tmp')


    Ires   = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = matrix(Ires$coefs,dim(coefs))

    write.table(ncoefs,file='optcoefs.tmp',col.names=FALSE,row.names=FALSE)
    write.table(ncoefs,file='curcoefs.tmp',col.names=FALSE,row.names=FALSE)

    if(out.meth=="optim")
    {
      if( is.null(control.out$meth) ){ control.out$meth = "BFGS" }
    	res = optim(pars[active],ProfileErr,allpars=pars,times=times,data=data,coef=coefs,lik=lik,proc=proc,active=active,
        hessian=T,in.meth=in.meth,control.in=control.in,control=control.out,gr=ProfileDP,method=control.out$meth)
    	npar = res$par
   	}

    else if(out.meth=="nlminb")
    {
        res = nlminb(pars[active],ProfileErr,allpars=pars,times=times,data=data,coef=coefs,lik=lik,proc=proc,
    	   in.meth=in.meth,control.in=control.in,control=control.out,gr=ProfileDP,active=active)
        npar = res$par
    }

    else if(out.meth=="maxNR")
    {
        if(is.null(control.out$print.level)){control.out$print.level = 2}
        if(is.null(control.out$iterlim)){control.out$iterlim = 100}
    
        res = maxNR(ProfileErr,start=pars[active],allpars=pars,times=times,data=data,coef=coefs,lik=lik,proc=proc,
    	   in.meth=in.meth,control.in=control.in,active=active,sgn=-1,
         grad=ProfileDP,print.level=control.out$print.level,iterlim=control.out$iterlim)
        npar = res$estimate
     }
     else if(out.meth=="subplex")
     {
        res = subplex(pars[active],ProfileErr,control=control.out,hessian=FALSE,
         allpars=pars,times=times,data=data,coef=coefs,lik=lik,proc=proc,
    	   in.meth=in.meth,control.in=control.in,active=active)
        npar = res$par
     }
     else if(out.meth=="trust")
     {
        if(is.null(control.out$rinit)){ control.out$rinit = 1}
        if(is.null(control.out$rmax)){ control.out$rinit = 100}
        if(is.null(control.out$parscale)){ control.out$parscale = abs(pars[active]) }
        if(is.null(control.out$iterlim)){ control.out$iterlim = 100}
     
        res = trust(ProfileList,pars[active],rinit=control.out$rinit,
        rmax=control.out$rmax,parscale=control.out$parscale,iterlim=control.out$iterlim,
          allpars=pars,times=times,data=data,coef=coefs,lik=lik,proc=proc,
    	    in.meth=in.meth,control.in=control.in,active=active)
        npar = res$argument
     }
     else if(out.meth == "ProfileGN"){
      res=Profile.GausNewt(pars=pars,times=times,data=data,coefs=ncoefs,
		    lik=lik,proc=proc,in.meth=in.meth,control.in=control.in,
		    active=active,control=control.out)
      npar = res$pars[active]

      ncoefs = res$in.res$coefs
      g = res$in.res$df
      resid = res$in.res$f
    }
    else if(out.meth == "nls"){
      res = nls(~ProfileSSE(pars,allpars,times,data,coefs,lik,proc,in.meth,control.in,active),
        data = list(allpars=pars,times=times,data=data,coefs=ncoefs,lik=lik,proc=proc,
        in.meth=in.meth,control.in=control.in,active=active),
        start = list(pars=pars[active]),trace=control.out$trace)
      npar = res$m$getPars()

      g = res$m$gradient()
      resid = res$m$resid()
    }
     
     if(file.exists('curcoefs.tmp')){
      	 ncoefs = as.matrix(read.table('curcoefs.tmp'))
   	 }
   	 else{ ncoefs = c() }
     if(file.exists('counter.tmp')){
      counter = as.matrix(read.table('counter.tmp'))
     } 	 
     else{  counter = c() }
    
    if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
    if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
    if(file.exists('counter.tmp')){file.remove('counter.tmp')}
        
    pars[active] = npar
    if(!is.null(proc$more$names)){ colnames(ncoefs) = proc$more$names }
    
    return(list(pars=pars,coefs=ncoefs,res=res,counter=counter))
}

##############################################################################

FitMatchOpt <- function(coefs,which,pars,proc,meth='nlminb',control=NULL)
{

   if(meth=="optim"){
        res = optim(coefs[,which],FitMatchErr,gr=FitMatchDC,hessian=T,
                    control=control,allcoefs=coefs,which=which,pars=pars,
                    proc=proc,method="BFGS")
        ncoefs = matrix(res$par,length(which),length(res$par)/length(which))
    }
    else if(meth=="nlminb"){
        res = nlminb(start=coefs[,which],objective=FitMatchErr,gradient=FitMatchDC,
                     hessian=FitMatchDC2,control=control,allcoefs=coefs,which=which,
                     pars=pars,proc=proc)
        ncoefs = matrix(res$par,length(which),length(res$par)/length(which))
    }
    else if(meth=="maxNR"){
        res = maxNR(FitMatchErr,coefs[,which],allcoefs=coefs,which=which,pars=pars,
                    proc=proc,sgn=-1,grad=FitMatchDC,hess=FitMatchDC2,
                    print.level=control$print.level,iterlim = control$iterlim)
        ncoefs = matrix(res$estimate,length(which),length(res$estimate)/length(which))
    }
    else if(meth=='trust'){
        if(is.null(control$rinit))   control$rinit = 1
        if(is.null(control$rmax))    control$rinit = 100
        if(is.null(control$parscale))control$parscale = length(length(which)*nrow(coefs))
        if(is.null(control$iterlim)) control$iterlim = 100
        res = trust(FitMatchList,coefs[,which],rinit=control$rinit,rmax=control$rmax,
                    parscale=control$parscale,iterlim=control$iterlim,
                    allcoefs=coefs,which=which,pars=pars,proc=proc)
        ncoefs = matrix(res$argument,length(which),length(res$estimate)/length(which))
    }
    coefs[,which] = ncoefs

    return(list(coefs=coefs,res=res))
}

##############################################################################

ParsMatchOpt <- function(pars,coefs,proc,active=1:length(pars),meth='nlminb',control=NULL)
{

   if(meth=="optim"){
        res = optim(pars[active],ParsMatchErr,gr=ParsMatchDP,hessian=T,
            control=control,coefs=coefs,proc=proc,active=active,allpars=pars,method="BFGS")
        npars = res$par
    }
    else if(meth=="nlminb"){
        res = nlminb(pars[active],ParsMatchErr,gradient=ParsMatchDP,
            control=control,coefs=coefs,proc=proc,active=active,allpars=pars)
        npars = res$par
    }
    else if(meth=="maxNR"){
        res = maxNR(ParsMatchErr,pars[active],coefs=coefs,proc=proc,active=active,
                    allpars=pars,sgn=-1,grad=ParsMatchDP,print.level=control$print.level,
                    iterlim = control$iterlim)
        npars = res$estimate
    }
    else if(meth=="trust")
    {
        if(is.null(control$rinit)){ control$rinit = 1}
        if(is.null(control$rmax)){ control$rinit = 100}
        if(is.null(control$parscale)){ control$parscale = length(active)}
        if(is.null(control$iterlim)){ control$iterlim = 100}

        res = trust(ParsMatchList,pars[active],rinit=control$rinit,rmax=control$rmax,
                    parscale=control$parscale,iterlim=control$iterlim,
                    coefs=coefs,proc=proc,active=active,allpars=pars)
        npars = res$argument
    }
    
    
    pars[active] = npars

    return(list(pars=pars,res=res))
}
