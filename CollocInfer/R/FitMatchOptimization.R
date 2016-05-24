######
# This file contains 'FitMatch' functions. These allow the estimates for some
# states to be updated while holding others fixed so that the whole system 
# minimizes the objective 'proc'; usually this measures departures from 
# a differential equation. 
#
# This will typically be used to obtain starting values for the inner optimization
# when some state variables are measured and can be smoothed, but other state variables
# are unmeasured. They are then estimated using FitMatchOpt.
#
# The file includes functions
#
# FitMatchOpt -- this calls various optimization routines
# FitMatch* -- gives the optimization criteria and its derivatives.  

FitMatchOpt <- function(coefs,which,pars,proc,meth='nlminb',control=list())
{
    check.lik.proc.data.coefs(lik=NULL,proc,data=NULL,times=NULL,coefs=coefs)

   if(meth=="optim"){
      if( is.null(control$trace) ){ control$trace = 6 }
      if( is.null(control$maxit) ){ control$maxit = 1000 }
      if( is.null(control$reltol) ){ control$reltol = 1e-12}
      if( is.null(control$meth) ){ control$meth = "BFGS" }
        res = optim(coefs[,which],FitMatchErr,gr=FitMatchDC,hessian=T,
            control=control,allcoefs=coefs,which=which,pars=pars,proc=proc,method="BFGS")

        ncoefs = matrix(res$par,length(which),length(res$par)/length(which))
    }
    else if(meth=="nlminb"){
      if( is.null(control$trace) ){ control$trace = 10 }
      if( is.null(control$eval.max) ){ control$eval.max = 2000 }
      if( is.null(control$iter.max) ){ control$iter.max = 1000 }
      if( is.null(control$rel.tol) ){ control$rel.tol = 1e-12 }       
        res = nlminb(start=coefs[,which],objective=FitMatchErr,gradient=FitMatchDC,hessian=FitMatchDC2,
            control=control,allcoefs=coefs,which=which,pars=pars,proc=proc)

        ncoefs = matrix(res$par,length(which),length(res$par)/length(which))
    }
    else if(meth=="maxNR"){
        if(is.null(control$print.level)){control$print.level = 2}
        if(is.null(control$iterlim)){control$iterlim = 1000}
        if(is.null(control$reltol)){control$reltol = 1e-12}
        res = maxLik::maxNR(FitMatchErr,coefs[,which],allcoefs=coefs,which=which,pars=pars,proc=proc,sgn=-1,
            grad=FitMatchDC,hess=FitMatchDC2,print.level=control$print.level,
            iterlim = control$iterlim)

        ncoefs = matrix(res$estimate,length(which),length(res$estimate)/length(which))
    }
    else if(meth=='trust'){
        if(is.null(control$rinit)){ control$rinit = 1}
        if(is.null(control$rmax)){ control$rmax = 100}
        if(is.null(control$parscale)){ control$parscale = length(length(which)*nrow(coefs))}
        if(is.null(control$iterlim)){ control$iterlim = 100}

        res = trust::trust(FitMatchList,coefs[,which],rinit=control$rinit,rmax=control$rmax,parscale=control$parscale,iterlim=control$iterlim,
          allcoefs=coefs,which=which,pars=pars,proc=proc)

        ncoefs = matrix(res$argument,length(which),length(res$estimate)/length(which))
    }
    else{
      stop('Unknown optimizer specified')
    }
    
    coefs[,which] = ncoefs

    return(list(coefs=coefs,res=res))
}


################################################

FitMatchErr = function(coefs,allcoefs,which,pars,proc,sgn=1)       # Matching unmeasured components
{                                                                  # which represents columns of
    allcoefs[,which] = matrix(coefs,nrow(allcoefs),length(which))  # allcoefs to be optimized

    g = proc$fn(allcoefs,proc$bvals,pars,proc$more)

    return(sgn*g)
}

################################################

FitMatchDC = function(coefs,allcoefs,which,pars,proc,sgn=1)    
{
    # Derivatives for matching unmeasured components 
   
    allcoefs[,which] = matrix(coefs,nrow(allcoefs),length(which))
 
    g = matrix( proc$dfdc(allcoefs,proc$bvals,pars,proc$more), dim(allcoefs) )

    g = as.vector( g[,which] )

    return(sgn*g)
}

################################################

FitMatchDC2 = function(coefs,allcoefs,which,pars,proc,sgn=1)    # Hessian for matching unmeasured components
{
    allcoefs[,which] = matrix(coefs,nrow(allcoefs),length(which))
    
    ind = array(1:length(allcoefs),dim(allcoefs))
    ind = as.vector(ind[,which])

    H = proc$d2fdc2(allcoefs,proc$bvals,pars,proc$more)

    return(as.matrix(sgn*H[ind,ind]))
}

################################################

FitMatchList = function(coefs,allcoefs,which,pars,proc,sgn=1)  # Inner Objective
{
  value = FitMatchErr(coefs,allcoefs,which,pars,proc,sgn=1)
  gradient = FitMatchDC(coefs,allcoefs,which,pars,proc,sgn=1)
  hessian = FitMatchDC2(coefs,allcoefs,which,pars,proc,sgn=1)

  return(list(value=value,gradient=gradient,hessian=hessian))

}
