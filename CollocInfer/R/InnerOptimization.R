#####
# This file contains functions defining the inner optimization of the Profile
# estimation regime. It provides functions
#
# - inneropt: which calls various optimization routines 
# - SplineEst.NewtRaph: a quick and dirty Newton-Raphson scheme
# - SplineCoefs*: functions that define the inner objective function and its 
# derivatives. 
#####


inneropt <- function(data,times,pars,coefs,lik,proc,in.meth='nlminb',
                     control.in=list())
{              
#  The multivariate data in argment data are observations of one or more
#  functional variables or processes at time points in argument time.

#  The function optimizes the coefficients of the basis expansions for the
#  representations of the variables or processes in what is called the
#  inner optimization loop.  The parameter values are held fixed.

#  The fitting criterion for both the data and the differential equation
#  are user-defined in arguments lik and proc, respectively.

#  The optimization method for the inner loop is defined by argument in.meth.

#  The function returns the optimized coefficients and the residual values.

  check.lik.proc.data.coefs(lik,proc,data,times,coefs)

  if(in.meth=="SplineEst"){
    #  optimization by Newton-Raphson using function SplineEst.NewtRaph
      if(is.null(control.in$reltol)){ control.in$reltol = 1e-12 }
      if(is.null(control.in$maxit)){ control.in$maxit = 1000}
      if(is.null(control.in$maxtry)){ control.in$maxtry = 10 }
      if(is.null(control.in$trace)){ control.in$trace = 0 }
        res = SplineEst.NewtRaph(coefs,times,data,lik,proc,pars,control.in)
        ncoefs = matrix(res$coefs,ncol(lik$bvals),length(res$coefs)/ncol(lik$bvals))
    }
    else if(in.meth=="optim"){
        #  optimization using function optim which in turn calls function SplineCoefsErr
      if( is.null(control.in$trace) ){ control.in$trace = 0 }
      if( is.null(control.in$maxit) ){ control.in$maxit = 1000 }
      if( is.null(control.in$reltol) ){ control.in$reltol = 1e-12}
      if( is.null(control.in$meth) ){ control.in$meth = "BFGS" }
      if( is.null(control.in$reportHessian)){ control.in$reportHessian = TRUE }
      
      imeth = control.in$meth
      control.in$meth = NULL
      
        res = optim(coefs,SplineCoefsErr,gr=SplineCoefsDC,hessian=control.in$reportHessian,
            control=control.in,times=times,data=data,lik=lik,proc=proc,pars=pars,method=imeth)

        ncoefs = matrix(res$par,ncol(lik$bvals),length(res$par)/ncol(lik$bvals))
    }    
    else if(in.meth=="nlminb"){
    #  optimization using function nlminb
      if( is.null(control.in$trace) ){ control.in$trace = 0 }
      if( is.null(control.in$eval.max) ){ control.in$eval.max = 2000 }
      if( is.null(control.in$iter.max) ){ control.in$iter.max = 1000 }
      if( is.null(control.in$rel.tol) ){ control.in$rel.tol = 1e-12 }  
      if( is.null(control.in$useHessian) ){ Hessian = SplineCoefsDC2 }
      else{ Hessian = NULL } 
          
        res = nlminb(coefs,SplineCoefsErr,gradient=SplineCoefsDC,hessian=Hessian,
            control=control.in,times=times,data=data,lik=lik,proc=proc,pars=pars)

        ncoefs = matrix(res$par,ncol(lik$bvals),length(res$par)/ncol(lik$bvals))
    }
    else if(in.meth=="maxNR"){
    #  optimizaton using function maxNR
        if(is.null(control.in$print.level)){control.in$print.level = 1}
        if(is.null(control.in$iterlim)){control.in$iterlim = 1000}
        if(is.null(control.in$reltol)){control.in$reltol = 1e-12}
        if( is.null(control.in$useHessian) ){ Hessian = SplineCoefsDC2 }
        else{ Hessian = NULL } 
    
        res = maxLik::maxNR(SplineCoefsErr,coefs,times=times,data=data,lik=lik,proc=proc,pars=pars,sgn=-1,
            grad=SplineCoefsDC,hess=Hessian,print.level=control.in$print.level,
            iterlim = control.in$iterlim)

        ncoefs = matrix(res$estimate,ncol(lik$bvals),length(res$estimate)/ncol(lik$bvals))
    }
    else if(in.meth=='trust'){
    #  optimization using function trust
        if(is.null(control.in$rinit)){ control.in$rinit = 1}
        if(is.null(control.in$rmax)){ control.in$rmax = 100}
        if(is.null(control.in$parscale)){ control.in$parscale = abs(coefs)}
        if(is.null(control.in$iterlim)){ control.in$iterlim = 100}
        
        res = trust::trust(SplineCoefsList,coefs,rinit=control.in$rinit,
                    rmax=control.in$rmax,parscale=control.in$parscale,
                    iterlim=control.in$iterlim,
                    times=time,data=data,lik=lik,proc=proc,pars=pars,sgn=1)
    
        ncoefs = matrix(res$argument,dim(coefs))
    }
    else{
      stop('Unknown optimizer specified')
    }

    if(!is.null(proc$more$names)){ colnames(ncoefs) = proc$more$names }
  
    return(list(coefs=ncoefs,res=res))
}

################  Newton-Raphson routine   #######################

SplineEst.NewtRaph = function(coefs,times,data,lik,proc,pars,    
    control=list(reltol=1e-12,maxit=1000,maxtry=10,trace=0))    
{
    if(is.null(control)){ control = list() }
    if(is.null(control$reltol)){ control$reltol = 1e-12 }
    if(is.null(control$maxit)){ control$maxit = 1000}
    if(is.null(control$maxtry)){ control$maxtry = 10 }
    if(is.null(control$trace)){ control$trace = 0 }

  check.lik.proc.data.coefs(lik,proc,data,times,coefs)
  
    f0 = SplineCoefsErr(coefs,times,data,lik,proc,pars)
    g = SplineCoefsDC(coefs,times,data,lik,proc,pars)
    H = SplineCoefsDC2sparse(coefs,times,data,lik,proc,pars)

    eigs = eigen(H)$values
    gam = -2*min(eigs[!(eigs==0)])
    if(gam <0 ){ gam = -gam/4 }
    eye = diag(rep(1,ncol(H)))

 # gam = 0

 #   if(is.matrix(H)){ DC = -ginv(H+gam*eye)%*%g  }
 #   else{ DC = -as.matrix(solve(H+gam*eye,g)) }
 
    gradnorm1 = 1
    fundif = 1
    iter = 0
    f1 = f0

    coefs0 = as.vector(coefs)

    while( gradnorm1 > control$reltol & fundif > 0 & iter < control$maxit){
  
        iter = iter + 1    
        if( is.matrix(H) ){ DC = -ginv(H+gam*eye)%*%g }
        else{ DC = -as.matrix(solve(H+gam*eye,g)) }
   
        ntry = 0
        
        coefs1 = coefs0 
        
 #       if(control$trace >0){ print(c(,f0,mean(abs(DC)),0)) }

        while( (f1 >= f0) & (t(DC)%*%DC > control$reltol) & (ntry < control$maxtry) ){
        

            coefs1 = coefs0 + DC
            f1 = SplineCoefsErr(coefs1,times,data,lik,proc,pars)
 #         DC = DC/2            
            gam = gam*2
            if( is.matrix(H) ){ DC = -ginv(H+gam*eye)%*%g }
            else{ DC = -as.matrix(solve(H+gam*eye,g)) }           
            ntry = ntry + 1
  #          print(c(f1,f0, t(DC)%*%DC,control$reltonl,ntry,control$maxtry))
  #          print((f1 >= f0) & (t(DC)%*%DC > control$reltol) & (ntry < control$maxtry))
            
        }

        gam = gam/4
        coefs0 = coefs1
        
        g = SplineCoefsDC(coefs0,times,data,lik,proc,pars)
        H = SplineCoefsDC2sparse(coefs0,times,data,lik,proc,pars)    

        gradnorm1 = mean(abs(DC))
        fundif = (f0-f1)/abs(f0)
  
        f0 = f1      
        if(control$trace>1){ print(c(iter,ntry,f0,gradnorm1,fundif)) }

    }
    if(control$trace>0){ print(c(f0,gradnorm1,iter,gam)) }
     
    return(list(coefs=coefs0,g=g,value=f0,H=H))
}


###############################################################

SplineCoefsList = function(coefs,times,data,lik,proc,pars,sgn=1)  # Inner Objective
{                              
  #  called in function inneropt by function trust
  value = SplineCoefsErr(coefs,times,data,lik,proc,pars,sgn)
  gradient = SplineCoefsDC(coefs,times,data,lik,proc,pars,sgn)
  hessian = SplineCoefsDC2(coefs,times,data,lik,proc,pars,sgn)

  return(list(value=value,gradient=gradient,hessian=hessian))
}

###############################################################

SplineCoefsErr = function(coefs,times,data,lik,proc,pars,sgn=1)  # Inner Objective
{
  #  evaluates the objective function for inner optimization
  #  called in function inneropt by functions optim, nlminb, maxNR

  coefs2 = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
  devals = as.matrix(lik$bvals %*% coefs2)

  colnames(devals) = proc$more$names  
    
  f = sum(lik$fn(data,times,devals,pars,lik$more)) + 
          proc$fn(coefs2,proc$bvals,pars,proc$more)

if(!is.null(proc$report)){ print(f) }    
  return(sgn*f)
}

###############################################################

SplineCoefsDC = function(coefs,times,data,lik,proc,pars,sgn=1)  
{
    #  Inner gradient with respect to coefficients
    coefs2 = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals%*%coefs2)
    colnames(devals) = proc$more$names

    g = as.matrix(t(lik$bvals)%*%lik$dfdx(data,times,devals,pars,lik$more)) + 
        proc$dfdc(coefs2,proc$bvals,pars,proc$more)

    g = as.vector(g)
    return(sgn*g)
}

###############################################################

SplineCoefsDP = function(coefs,times,data,lik,proc,pars,sgn=1)  # Inner gradient
{
    #  Inner gradient with respect to parameters
    coefs2 = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals%*%coefs2)
    colnames(devals) = proc$more$names
    
    g = apply(lik$dfdp(data,times,devals,pars,lik$more),2,sum)    + 
        proc$dfdp(coefs2,proc$bvals,pars,proc$more)

    g = as.vector(g)
    return(sgn*g)
}

###############################################################

SplineCoefsDC2sparse = function(coefs,times,data,lik,proc,pars,sgn=1) 
{
  #  Inner Hessian with respect to coefficients
    coefs2 = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals%*%coefs2)
    colnames(devals) = proc$more$names
    
    d2lik = lik$d2fdx2(data,times,devals,pars,lik$more)

    H = list(len=ncol(lik$bvals))
    for(i in 1:dim(d2lik)[2]){
      H[[i]] = list(len=ncol(devals))
        for(j in 1:dim(d2lik)[3]){
            H[[i]][[j]] = t(lik$bvals)%*%diag(d2lik[,i,j])%*%lik$bvals
        }
    }

    H = blocks2mat(H) 

    H = H + proc$d2fdc2(coefs2,proc$bvals,pars,proc$more)

    return(sgn*H)
}

SplineCoefsDC2 = function(coefs,times,data,lik,proc,pars,sgn=1)
{
  result = as.matrix(SplineCoefsDC2sparse(coefs,times,data,lik,proc,pars,sgn))
  return(result)
}

###############################################################

SplineCoefsDCDP = function(coefs,times,data,lik,proc,pars,sgn=1) 
{
    # Inner cross derivatives with respect to coefficients and parameters
    coefs2 = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals%*%coefs2)
    colnames(devals) = proc$more$names
    
    d2lik = lik$d2fdxdp(data,times,devals,pars,lik$more)

    H = c()
    for(i in 1:length(pars)){
        H = cbind(H,as.vector(as.matrix(t(lik$bvals)%*%d2lik[,,i])))
    }

    H = H + proc$d2fdcdp(coefs2,proc$bvals,pars,proc$more)
    
    return( as.matrix(sgn*H ))
}

