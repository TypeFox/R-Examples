#####
# This file contains functions defining the inner optimization of the Profile
# estimation regime. It provides functions
#
# - outeropt: which calls various optimization routines 
# - Profile.GausNewt: a quick and dirty Gauss-Newton scheme
# - ProfileErr*: functions that define the inner objective function and its 
#                derivatives. 
# - Profile.covariance: functions to estimate sample covariance of profiled 
# parameters. 
#####

outeropt <- function(data, times, pars, coefs, lik, proc,
                     in.meth='nlminb',  out.meth='nlminb',
                     control.in=list(), control.out=list(),
                     active=1:length(pars))
{
#  The multivariate data in argment data are observations of one or more
#  functional variables or processes at time points in argument times.

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

#  The function returns the optimized parameter values, a functional data object  
#  definingthe smoothing functions, the likelihood and process functions, and 
#  the residual values.

    check.lik.proc.data.coefs(lik,proc,data,times,coefs)
    
    #  clear any old temporary files 

    if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
    if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
    if(file.exists('counter.tmp')) {file.remove('counter.tmp')}

    #  initial inner optimization to get coefficients
                                                                                                  
    Ires   = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = matrix(Ires$coefs,dim(coefs))

    #  write coefficients to temporary files

    write.table(ncoefs,file='optcoefs.tmp',col.names=FALSE,row.names=FALSE)
    write.table(ncoefs,file='curcoefs.tmp',col.names=FALSE,row.names=FALSE)

    #  Select optimization method and optimize with respect to parameters
    
    ##########  method = "optim"  ###########

    if(out.meth=="optim"){
      if( is.null(control.out$trace) ) { control.out$trace   = 6     }
      if( is.null(control.out$maxit) ) { control.out$maxit  = 100    }
      if( is.null(control.out$reltol) ){ control.out$reltol = 1e-8   }
      if( is.null(control.out$meth) )  { control.out$meth   = "BFGS" }
      
      ometh = control.out$meth
      control.out$meth = NULL
      
    	res = optim(pars[active], ProfileErr, allpars=pars, times=times,
                  data=data, coef=coefs, lik=lik, proc=proc, 
                  active=active, hessian=T, in.meth=in.meth, 
                  control.in=control.in, control=control.out,
                  gr=ProfileDP, method=ometh)
    	npar = res$par
   	}

    ##########  method = "nlminb"  ###########

    else if(out.meth=="nlminb")
    {
      if( is.null(control.out$trace) )   { control.out$trace    = 10   }
      if( is.null(control.out$eval.max) ){ control.out$eval.max = 200  }
      if( is.null(control.out$iter.max) ){ control.out$iter.max = 100  }
      if( is.null(control.out$rel.tol) ) { control.out$rel.tol  = 1e-8 }       
        
        res = nlminb(pars[active], ProfileErr, allpars=pars, times=times,
                     data=data, coef=coefs, lik=lik, proc=proc,
    	               in.meth=in.meth, 
                     control.in=control.in, control=control.out,
                     gradient=ProfileDP, active=active)
        npar = res$par
    }

    ##########  method = "maxNR"  ###########
    
    else if(out.meth=="maxNR")
    {
        if(is.null(control.out$print.level)){control.out$print.level = 2   }
        if(is.null(control.out$iterlim))    {control.out$iterlim     = 100 }
        if(is.null(control.out$reltol))     {control.out$reltol      = 1e-8}
    
        res = maxLik::maxNR(ProfileErr1, start=pars[active], allpars=pars, times=times,
                    data=data, coef=coefs, lik=lik, proc=proc,
    	              in.meth=in.meth, control.in=control.in, sgn=-1, 
                    active1=active, grad=ProfileDP1, 
                    print.level=control.out$print.level,
                    iterlim=control.out$iterlim)
        npar = res$estimate
     }

    ##########  method = "subplex"  ###########
    
     else if(out.meth=="subplex")
     {
        if( is.null(control.out$maxit)) { control.out$maxit  = 100  }
        if( is.null(control.out$reltol)){ control.out$reltol = 1e-8 }
        
        res = subplex::subplex(pars[active], ProfileErr, control=control.out,
                      hessian=FALSE, allpars=pars, times=times, data=data,
                      coef=coefs, lik=lik, proc=proc,
    	                in.meth=in.meth, control.in=control.in, active=active)
        npar = res$par
     }
     
    ##########  method = "trust"  ###########
    
#     else if(out.meth=="trust")
#     {
#        if(is.null(control.out$rinit))   { control.out$rinit    = 1}
#        if(is.null(control.out$rmax))    { control.out$rmax     = 100}
#        if(is.null(control.out$iterlim)) { control.out$itelim   = 100}
#        if(is.null(control.out$parscale)){ 
#            control.out$parscale = abs(pars[active]) 
#        }
#     
#        res = trust(ProfileList, pars[active], 
#                    rinit=control.out$rinit, rmax=control.out$rmax,
#                    parscale=control.out$parscale, iterlim=control.out$iterlim,
#                    allpars=pars,times=times, data=data, coef=coefs, 
#                    lik=lik, proc=proc,
#    	              in.meth=in.meth, control.in=control.in, active=active)
#        npar = res$argument
#     }
     
    ##########  method = "ProfileGN"  ###########
    
     else if(out.meth == "ProfileGN"){
      if(is.null(control.out$reltol)){ control.out$reltol = 1e-12 }
      if(is.null(control.out$maxit)) { control.out$maxit  = 50    }
      if(is.null(control.out$maxtry)){ control.out$maxtry = 15    }
      if(is.null(control.out$trace)) { control.out$trace  = 1     }
      res = Profile.GausNewt(pars=pars, times=times, data=data, coefs=ncoefs,
		               lik=lik, proc=proc, in.meth=in.meth, control.in=control.in,
                   active=active, control=control.out)
      npar = res$pars[active]

      ncoefs = res$in.res$coefs
      g      = res$in.res$df
      resid  = res$in.res$f
    }
    
    ##########  method = "nls"  ###########
    
    else if(out.meth == "nls"){
      if(is.null(control.out$trace))    {control.out$trace     = TRUE}
      if(is.null(control.out$maxiter))  {control.out$maxiter   = 100 }
      if(is.null(control.out$tol))      {control.out$tol       = 1e-8}
      if(is.null(control.out$printEval)){control.out$printEval = TRUE}
      if(is.null(control.out$warnOnly)) {control.out$warnOnly  = TRUE}
      
      res = nls(~ProfileSSE(pars, allpars, times, data, coefs, lik, proc,
                            in.meth, control.in, active),
                data = list(allpars=pars, times=times, data=data, coefs=ncoefs,
                            lik=lik, proc=proc, in.meth=in.meth,
                            control.in=control.in,active=active),
                start = list(pars=pars[active]), 
                trace = control.out$trace,
                control=control.out)
      npar  = res$m$getPars()
      g     = res$m$gradient()
      resid = res$m$resid()
    }
    else {
      stop("Unrecognized optimization method.")
    }
     
     if(file.exists('curcoefs.tmp')){
      	 ncoefs = as.matrix(read.table('curcoefs.tmp'))
   	 }
   	 else{ ncoefs = c() }
     if(file.exists('counter.tmp')){
      counter = as.matrix(read.table('counter.tmp'))
     } 	 
     else{  counter = c() }
    
    #  remove temporary files
    
    if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
    if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
    if(file.exists('counter.tmp')) {file.remove('counter.tmp') }
        
    #  return results
        
    pars[active] = npar
    if(!is.null(proc$more$names)){ colnames(ncoefs) = proc$more$names }
    
    return(list(pars=pars, coefs=ncoefs, res=res, counter=counter))
}

########################################################################


Profile.GausNewt = function(pars,times,data,coefs,lik,proc,   
    in.meth='nlminb',control.in=NULL,active=1:length(pars),   
    control=list(reltol=1e-6,maxit=50,maxtry=15,trace=1)) 
{

    # Gauss-Newton routine for squared error outer objective
    
    #  set defaults for control (control.out in call)
    if(is.null(control))       { control = list()       }
    if(is.null(control$reltol)){ control$reltol = 1e-12 }
    if(is.null(control$maxit)) { control$maxit = 50     }
    if(is.null(control$maxtry)){ control$maxtry = 15    }
    if(is.null(control$trace)) { control$trace = 1      }

    #  check lik and proc objects
    
    check.lik.proc.data.coefs(lik,proc,data,times,coefs)

    #  Initial inner optimization
    
    res0 = ProfileSSE.AllPar(pars, times, data, coefs, 
              lik, proc, in.meth, control.in, use.nls=FALSE)

     #  allows some parameters to be fixed

    ind = names(pars)
    if(is.null(ind)){ ind = 1:length(pars) }                
    ind = ind[active]                                      

    pars0 = pars
    pars1 = pars
    F0 = sum(res0$value^2)
    F1 = F0
    iter = 0
    fundif = 1


    Jacobian = res0$gradient[,active]
    residual = res0$value

    gradnorm0 = mean(abs(crossprod(Jacobian, residual)))
    gradnorm1 = 1
    dcdp = c()

    eigs = eigen(crossprod(Jacobian,Jacobian))$values
    gam = -2*min(eigs[!(eigs==0)])
    if(gam <0 ){ gam = -gam/4 }
    eye = diag(rep(1,ncol(Jacobian)))

#    gam = 0

    while( gradnorm1 > control$reltol & 
           fundif    > control$reltol & 
           iter      < control$maxit ) {
  
        iter = iter + 1

 #       Dpars = lsfit(Jacobian, residual, intercept=FALSE)
        Dpars = solve(crossprod(Jacobian)+gam*eye, crossprod(Jacobian,residual))
    
        ntry = 0
        
        while(F1 >= F0 & 
              t(Dpars)%*%Dpars > control$reltol & 
              ntry < control$maxtry){
              
              
            pars1[active] = pars0[active] - Dpars
            tres = ProfileSSE.AllPar(pars1, times, data, res0$coefs, lik, proc,
                                     in.meth, control.in, res0$dcdp, pars0,
                                     use.nls=FALSE)
            F1    = sum(tres$value^2)
#           Dpars = Dpars/2
            gam = gam*2
            Dpars = solve(crossprod(Jacobian)+gam*eye, 
                          crossprod(Jacobian,residual))
            ntry  = ntry + 1
        }

        gam = gam/4
        res0      = tres
        dcdp      = tres$dcdp

        gradnorm1 = abs(mean(t(res0$gradient[,active])%*%res0$value))
        fundif    = (F0-F1)/abs(F0)

        pars0     = pars1
        res0      = tres
        F0        = F1
        gradnorm0 = gradnorm1
        
        if(control$trace>0){ print(c(iter,ntry,F0,pars0)) }
    }

    newpars = pars0
     
    return(list(pars=newpars,in.res=res0,value=F0))
}

##################################################################

ProfileErr.AllPar = function(pars, times, data, coefs, lik, proc,
                             in.meth = "house", control.in=NULL, sgn=1)
{                                                          
#    coefs = as.vector(coefs)  # First run the inner optimization
    
    if(file.exists('curcoefs.tmp')) {                      
      altcoefs = as.matrix(read.table('curcoefs.tmp'))               
      if( !(length(altcoefs)==length(coefs)) ){
        stop(paste('Variables in curcoefs.tmp do not conform;',
                   ' file exists from previous experiments?'))
       }
    }
    else{
      altcoefs = coefs
    }
    
    if(file.exists('counter.tmp')){                       
      counter = read.table('counter.tmp')                 
      niter = counter[nrow(counter),1]
    }
    else{
      counter = matrix(c(1,0,pars),1,length(pars)+2)
      niter = 0
    }
    
           
    altdevals = as.matrix(lik$bvals%*%matrix(altcoefs, ncol(lik$bvals),
                          length(altcoefs)/ncol(lik$bvals)))

    colnames(altdevals) = proc$more$names
    
    f1 = SplineCoefsErr(coefs,times,data,lik,proc,pars)
    f2 = SplineCoefsErr(altcoefs,times,data,lik,proc,pars)    

    if(f2 < f1){ coefs = altcoefs }

    Ires   = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = Ires$coefs

    devals = as.matrix(lik$bvals %*% ncoefs)
    colnames(devals) = proc$more$names
    
    f  = sum(lik$fn(data,times,devals,pars,lik$more))

    f2 = sum(lik$fn(data,times,altdevals,pars,lik$more))

    #    wout = 0
    if(f <= f2){
       write.table(ncoefs,file='curcoefs.tmp',col.names=FALSE,row.names=FALSE)  
    #      save(ncoefs,niter,pars,file='bestcoefs.Rdata')
    }                                                                   
    write.table(ncoefs,file='optcoefs.tmp',col.names=FALSE,row.names=FALSE)

    if(niter==0){ 
      counter[1,2] = f
      write.table(counter,file='counter.tmp',col.names=FALSE,row.names=FALSE) 
    }
 
    if(niter>=1){
      if(f < counter[niter,2]){
       counter = rbind(counter,c(niter+1,f,pars))
       write.table(counter,file='counter.tmp',col.names=FALSE,row.names=FALSE)
      }
    }

     if(!is.null(lik$report)){ print(f) }
    f = sgn*f
    
#    if(file.exists('ProfHistory.Rdata')){ load('ProfHistory.Rdata') }
#    else{ 
#      parmat = c()
#      fvec = c()
#      coefmat = c()    
#      wvec = c()
#    }
#    
#    parmat = cbind(parmat,pars)
#    coefmat = cbind(coefmat,as.vector(ncoefs))
#    fvec = c(fvec,f)
#    wvec = c(wvec,wout)
#    
#    save(parmat,coefmat,fvec,wvec,file='ProfHistory.Rdata')
    

#    print(c(f,pars))
    return(f)
}

##################################################################

ProfileDP.AllPar = function(pars, times, data, coefs, lik, proc,
                            in.meth=NULL, control.in=NULL, sgn=1, sumlik=1)
{

    if(file.exists('optcoefs.tmp')){
      altcoefs = as.matrix(read.table('optcoefs.tmp'))
      if( !(length(altcoefs)==length(coefs)) ){
        stop(paste('Variables in curcoefs.tmp do not conform;',
                   'file exists from previous experiments?'))
      }else{
        coefs = altcoefs
      }      
    }

    coefs  = as.matrix(coefs)
    devals = as.matrix(lik$bvals%*%coefs)
    colnames(devals) = proc$more$names

#    coefs = as.vector(coefs)    

    d2Hdc2  = SplineCoefsDC2sparse(coefs,times,data,lik,proc,pars)
    d2Hdcdp = SplineCoefsDCDP(coefs,times,data,lik,proc,pars)
  
    if(is.matrix(d2Hdc2)){ 
        dcdp = -ginv(d2Hdc2)%*%d2Hdcdp 
    } else{ 
        dcdp = -as.matrix(solve(d2Hdc2,d2Hdcdp)) 
    }
  
   
    if(sumlik){
        dlikdc = as.vector(t(lik$bvals) %*% 
                           lik$dfdx(data, times, devals, pars, lik$more))
        df = t(dcdp) %*% (dlikdc) + 
             apply(lik$dfdp(data, times, devals, pars, lik$more),2,sum)    
    
        df = as.vector(sgn*df)
        return(df)
    } else{
        dlikdx = lik$dfdx(data,times,devals,pars,lik$more)

        dlikdp = lik$dfdp(data,times,devals,pars,lik$more)

        dlikdc = c()
        for(i in 1:dim(dlikdx)[2]){
            dlikdc = cbind(dlikdc,as.matrix(diag(dlikdx[,i])%*%lik$bvals))            
        }       
        
        df = dlikdc%*%dcdp + dlikdp
        
        return(as.matrix(df))
    }
}

##################################################################

ProfileErr = function(pars,allpars,times,data,coefs,lik,proc,in.meth = "house",
        control.in=NULL,sgn=1,active=1:length(allpars))
{
#print(allpars)
#print(active)
    allpars[active] = pars
#print(allpars)

    f = ProfileErr.AllPar(allpars, times, data, coefs, lik, proc,
                          in.meth, control.in,sgn)
    return(f)
}


# The following fixes a conflict that arises when using maxNR. 
ProfileErr1 =    function(pars,allpars,times,data,coefs,lik,proc,
                          in.meth = "house",control.in=NULL,sgn=1,
                          active1=1:length(allpars))
{
  ProfileErr(pars,allpars,times,data,coefs,lik,proc,in.meth = "house",
        control.in=NULL,sgn=1,active=active1)
}

##################################################################

ProfileDP = function(pars, allpars, times, data, coefs, lik, proc,
                     in.meth = "house", control.in=NULL, sgn=1, sumlik=1,
                     active=1:length(allpars))
{
    allpars[active] = pars

    g = ProfileDP.AllPar(allpars, times, data, coefs, lik, proc,
                         in.meth, control.in, sgn, sumlik)


    if(sumlik){ 
        names(g) = names(allpars)
        g = g[active] 
    }
    else{ 
        colnames(g) = names(allpars)
        g = g[,active,drop=FALSE] 
    }
    
    return(g)        
}

##################################################################

ProfileDP1 = function(pars, allpars, times, data, coefs, lik, proc,
                     in.meth = "house", control.in=NULL, sgn=1, sumlik=1,
                     active1=1:length(allpars))
{                     
  ProfileDP(pars, allpars, times, data, coefs, lik, proc,
                     in.meth = "house", control.in=NULL, sgn=1, sumlik=1,
                     active=active1)                     
}                    
##################################################################

ProfileList = function(pars, allpars, times, data, coefs, lik, proc,
                       in.meth = "house", control.in=NULL, sgn=1, 
                       active=1:length(allpars))
{
  check.lik.proc.data.coefs(lik,proc,data,times,coefs)

  value = ProfileErr(pars=pars, allpars=allpars, times=times, data=data, 
                       coefs=coefs, lik=lik, proc=proc,
                       in.meth=in.meth, control.in=control.in, sgn=sgn, 
                       active=active)
  gradient = ProfileDP(pars=pars, allpars=allpars, times=times, data=data, 
                       coefs=coefs, lik=lik, proc=proc,
                       in.meth=in.meth, control.in=control.in, sgn=sgn, 
                       active=active)

  return(list(value=value,gradient=gradient))
}

##################################################################

ProfileSSE = function(pars, allpars, times, data, coefs, lik, proc,
                      in.meth='nlminb', control.in=NULL,   
                      active=1:length(pars), dcdp=NULL, oldpars=NULL,
                      use.nls=TRUE, sgn=1)             
{
  # Squared Error outer criterion
  
  allpars[active] = pars

  f = ProfileSSE.AllPar(pars=allpars, times=times, data=data, coefs=coefs,
                        lik=lik, proc=proc, in.meth=in.meth,
                        control.in=control.in, dcdp=dcdp, oldpars=oldpars,
                        use.nls=use.nls, sgn=sgn)
  
    if(use.nls){
        attr(f,'gradient') = attr(f,'gradient')[,active]
    }
    else{
        f$gradient = f$gradient[,active,drop=FALSE]
        f$dcdp = f$dcdp[,active,drop=FALSE]
    }
    return(f)
}

##################################################################

ProfileSSE.AllPar = function(pars, times, data, coefs, lik, proc,
                             in.meth='nlminb', control.in=NULL,   
                             dcdp=NULL, oldpars=NULL, use.nls=TRUE, sgn=1)                            
{           
  # Squared Error outer criterion

    #    coefs = as.vector(coefs)  # First run the inner optimization
    f1 = SplineCoefsErr(coefs,times,data,lik,proc,pars)

    if(use.nls){  
        # If using NLS, we need to keep track                                      
        if(file.exists('curcoefs.tmp')){
          altcoefs = as.matrix(read.table('curcoefs.tmp'))
          if( !(length(altcoefs)==length(coefs)) ){
            stop(paste('Variables in curcoefs.tmp do not conform;',
                       'file exists from previous experiments?'))
          }
        } else {
          altcoefs = coefs
        }
        if(file.exists('counter.tmp')){                       
          counter = read.table('counter.tmp')                 
          niter = counter[nrow(counter),1]
        } else {
          counter = matrix(c(1,0,pars),1,length(pars)+2)
          niter = 0
        }

        f2 = SplineCoefsErr(altcoefs,times,data,lik,proc,pars)    
    
        if(f2 < f1){
          coefs = altcoefs
          f1 = f2
        }
        
        altdevals = as.matrix(lik$bvals %*% matrix(altcoefs,ncol(lik$bvals),
                              length(altcoefs)/ncol(lik$bvals)))
        colnames(altdevals) = proc$more$names
    }

    if( !is.null(dcdp) ){
        tcoefs = as.vector(coefs) + dcdp%*%(pars-oldpars);
        f2 = SplineCoefsErr(tcoefs,times,data,lik,proc,pars)
        if(f2 < f1){
          coefs = tcoefs
          f1 = f2
        }         
    }

    Ires = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = Ires$coefs
    
    devals = as.matrix(lik$bvals%*%ncoefs)
    colnames(devals) = proc$more$names


    weights = checkweights(lik$more$weights,lik$more$whichobs,data)
    f = as.vector( as.matrix(data - 
        lik$more$fn(times, devals, pars, lik$more$more))*sqrt(weights) )
    isnaf = is.na(f)
    f[isnaf] = 0

    # Now a gradient

    dlikdp = lik$more$dfdp(times,devals,pars,lik$more$more)
    dlikdp = matrix(dlikdp,dim(dlikdp)[1]*dim(dlikdp)[2],dim(dlikdp)[3])

    dlikdx = lik$more$dfdx(times,devals,pars,lik$more$more)

    dlikdc = c()
    for(i in 1:dim(dlikdx)[2]){
        tH = c()
        for(j in 1:dim(dlikdx)[3]){
            tH = cbind(tH,as.matrix(diag(dlikdx[,i,j])%*%lik$bvals))            
        }
        dlikdc = rbind(dlikdc,tH)
    }

    d2Hdc2  = SplineCoefsDC2sparse(ncoefs,times,data,lik,proc,pars)
    d2Hdcdp = SplineCoefsDCDP(ncoefs,times,data,lik,proc,pars)

    if(is.matrix(d2Hdc2)){ 
        dcdp = ginv(d2Hdc2) %*% d2Hdcdp 
    } else { 
        dcdp = as.matrix(solve(d2Hdc2,d2Hdcdp)) 
    }

    df = dlikdc%*%dcdp + dlikdp
    df[isnaf,] = 0
    colnames(df) = proc$more$parnames

    if(!is.null(lik$report)){ print(f) }

    f = sgn*f
    df = sgn*df

    if(use.nls){
        tf = sum(lik$fn(data,times,devals,pars,lik$more))
        tf2 = sum(lik$fn(data,times,altdevals,pars,lik$more))
                                                                                                      
        if(tf <= tf2){
            write.table(ncoefs,file='curcoefs.tmp',
                        row.names=FALSE,col.names=FALSE)
            niter = counter[nrow(counter),1]
        }
        
        if(niter==0){ 
          counter[1,2] = tf
          write.table(counter,file='counter.tmp',
                      col.names=FALSE,row.names=FALSE) 
        }

     if(niter > 1){
      if(tf < counter[niter,2]){
       counter = rbind(counter,c(niter+1,tf,pars))
       write.table(counter,file='counter.tmp',col.names=FALSE,row.names=FALSE)
      }
    }
        
        write.table(ncoefs,file='optcoefs.tmp',col.names=FALSE,row.names=FALSE)       
        attr(f,'gradient') = df
        return(f)
    }
    else{
        return(list(value=f,gradient=df,coefs=ncoefs,dcdp=dcdp))
    }
}
