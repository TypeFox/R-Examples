######
# This file contains 'ParsMatch' functions. These allow initial 'gradient matching'
# estimates of parameters using fixed estimates for the state variables. 
#
# This will typically be used to obtain starting values for parameters when
# all state variables are measured and can be smoothed. 
#
# The file includes functions
#
# ParsMatchOpt -- this calls various optimization routines
# ParsMatch* -- gives the optimization criteria and its derivatives.  

ParsMatchOpt <- function(pars,coefs,proc,active=1:length(pars),meth='nlminb',control=list())
{
  check.lik.proc.data.coefs(lik=NULL,proc,data=NULL,times=NULL,coefs=coefs)

   if(meth=="optim"){
      if( is.null(control$trace) ){ control$trace = 6 }
      if( is.null(control$maxit) ){ control$maxit = 100 }
      if( is.null(control$reltol) ){ control$reltol = 1e-8}
      if( is.null(control$meth) ){ control$meth = "BFGS" }
        res = optim(pars[active],ParsMatchErr,gr=ParsMatchDP,hessian=T,
            control=control,coefs=coefs,proc=proc,active=active,allpars=pars,method="BFGS")

        npars = res$par
    }
    else if(meth=="nlminb"){
      if( is.null(control$trace) ){ control$trace = 10 }
      if( is.null(control$eval.max) ){ control$eval.max = 200 }
      if( is.null(control$iter.max) ){ control$iter.max = 100 }
      if( is.null(control$rel.tol) ){ control$rel.tol = 1e-8 }       
        res = nlminb(pars[active],ParsMatchErr,gradient=ParsMatchDP,
            control=control,coefs=coefs,proc=proc,active=active,allpars=pars)

        npars = res$par
    }
    else if(meth=="maxNR"){
        if(is.null(control$print.level)){control$print.level = 2}
        if(is.null(control$iterlim)){control$iterlim = 100}
        if(is.null(control$reltol)){control$reltol = 1e-8}
        res = maxLik::maxNR(ParsMatchErr,pars[active],coefs=coefs,proc=proc,activePar=active,allpars=pars,sgn=-1,
            grad=ParsMatchDP,print.level=control$print.level,
            iterlim = control$iterlim)

        npars = res$estimate
    }
    else if(meth=="trust")
    {
        if(is.null(control$rinit)){ control$rinit = 1}
        if(is.null(control$rmax)){ control$rmax = 100}
        if(is.null(control$parscale)){ control$parscale = abs(pars[active]) }
        if(is.null(control$iterlim)){ control$itelim = 100}

        res = trust::trust(ParsMatchList,pars[active],rinit=control$rinit,rmax=control$rmax,parscale=control$parscale,iterlim=control$iterlim,
          coefs=coefs,proc=proc,active=active,allpars=pars)
        npars = res$argument
     }
    else{
      stop('Unknown optimizer specified')
    }    
    
    pars[active] = npars

    return(list(pars=pars,res=res))
}


###################################################

ParsMatchErr = function(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)    # "Gradient Matching": fix the state and minimize proc
{
   allpars[active] = pars
   
   f = proc$fn(coefs,proc$bvals,allpars,proc$more)

   return(sgn*f)
}

###################################################

ParsMatchDP = function(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)   # Derivative for Gradient Matching
{
    allpars[active] = pars

   g = proc$dfdp(coefs,proc$bvals,allpars,proc$more)
   
   return(sgn*g[active])
}

###################################################

ParsMatchList = function(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)  
{
# A list of ParsMatch objective and gradient for the 'trust' optimization routine
  value = ParsMatchErr(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)
  gradient = ParsMatchDP(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)

  return(list(value=value,gradient=gradient))

}
