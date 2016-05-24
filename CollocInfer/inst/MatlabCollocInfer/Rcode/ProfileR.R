#  This file contains many of the service routines required by
#  CollocInfer for optimization given a user-defined 
#  lik and proc functions.

###############################################################

SplineCoefsList = function(coefs,times,data,lik,proc,pars,sgn=1)  
{
  value    = SplineCoefsErr(coefs,times,data,lik,proc,pars,sgn)
  gradient = SplineCoefsDC(coefs,times,data,lik,proc,pars,sgn)
  hessian  = SplineCoefsDC2(coefs,times,data,lik,proc,pars,sgn)
  return(list(value=value,gradient=gradient,hessian=hessian))
}

###############################################################

SplineCoefsErr = function(coefs,times,data,lik,proc,pars,sgn=1)  
{
  coefs = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
  devals = as.matrix(lik$bvals%*%coefs)

  colnames(devals) = proc$more$names  
    
  f = sum(lik$fn(data,times,devals,pars,lik$more)) + 
          proc$fn(coefs,proc$bvals,pars,proc$more)
  if(!is.null(proc$report)) print(f)    
  return(sgn*f)
}

###############################################################

SplineCoefsDC = function(coefs,times,data,lik,proc,pars,sgn=1)  
{
    #  Inner gradient with respect to coefficients
    coefs  = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals%*%coefs)
    colnames(devals) = proc$more$names

    g = as.matrix(t(lik$bvals)%*%lik$dfdx(data,times,devals,pars,lik$more)) + 
        proc$dfdc(coefs,proc$bvals,pars,proc$more)

    g = as.vector(g)
    return(sgn*g)
}

###############################################################

SplineCoefsDP = function(coefs,times,data,lik,proc,pars,sgn=1)  
{
    #  Inner gradient with respect to parameters
    coefs = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals%*%coefs)
    colnames(devals) = proc$more$names
    
    g = apply(lik$dfdp(data,times,devals,pars,lik$more),2,sum)    + 
        proc$dfdp(coefs,proc$bvals,pars,proc$more)

    g = as.vector(g)
    return(sgn*g)
}

###############################################################

SplineCoefsDC2 = function(coefs,times,data,lik,proc,pars,sgn=1) 
{
    #  Inner Hessian with respect to coefficients
    coefs = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals%*%coefs)
    colnames(devals) = proc$more$names
    d2lik = lik$d2fdx2(data,times,devals,pars,lik$more)
    H = list(len=ncol(lik$bvals))
    for(i in 1:dim(d2lik)[2]){
      H[[i]] = list(len=ncol(devals))
        for(j in 1:dim(d2lik)[3]){
            H[[i]][[j]] = as.matrix(t(lik$bvals)%*%diag(d2lik[,i,j])%*%lik$bvals)
        }
    }

    H = blocks2mat(H) + proc$d2fdc2(coefs,proc$bvals,pars,proc$more)
    

    return( as.matrix(sgn*H) )
}

###############################################################

SplineCoefsDCDP = function(coefs,times,data,lik,proc,pars,sgn=1) 
{
    # Inner cross derivatives with respect to coefficients and parameters
    coefs = matrix(coefs,ncol(lik$bvals),length(coefs)/ncol(lik$bvals))
    devals = as.matrix(lik$bvals%*%coefs)
    colnames(devals) = proc$more$names
    
    d2lik = lik$d2fdxdp(data,times,devals,pars,lik$more)

    H = c()
    for(i in 1:length(pars)){
        H = cbind(H,as.vector(as.matrix(t(lik$bvals)%*%d2lik[,,i])))
    }

    H = H + proc$d2fdcdp(coefs,proc$bvals,pars,proc$more)
    
    return( as.matrix(sgn*H ))
}

###############################################################

ParsMatchErr = function(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)    # "Gradient Matching": fix the state and minimize proc
{
   allpars[active] = pars
   
   f = proc$fn(coefs,proc$bvals,allpars,proc$more)

   return(sgn*f)
}

###############################################################

ParsMatchDP = function(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)   # Derivative for Gradient Matching
{
    allpars[active] = pars

   g = proc$dfdp(coefs,proc$bvals,allpars,proc$more)
   
   return(sgn*g[active])
}

###############################################################

ParsMatchList = function(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)  # Inner Objective
{
  value = ParsMatchErr(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)
  gradient = ParsMatchDP(pars,coefs,proc,active=1:length(pars),allpars,sgn=1)

  return(list(value=value,gradient=gradient))

}

###############################################################

FitMatchErr = function(coefs,allcoefs,which,pars,proc,sgn=1)       # Matching unmeasured components
{                                                                  # which represents columns of
    allcoefs[,which] = matrix(coefs,nrow(allcoefs),length(which))  # allcoefs to be optimized

    g = proc$fn(allcoefs,proc$bvals,pars,proc$more)

    return(sgn*g)
}

###############################################################

FitMatchDC = function(coefs,allcoefs,which,pars,proc,sgn=1)    # Derivatives for matching unmeasured components
{
    allcoefs[,which] = matrix(coefs,nrow(allcoefs),length(which))

    g = matrix( proc$dfdc(allcoefs,proc$bvals,pars,proc$more), dim(allcoefs) )

    g = as.vector( g[,which] )

    return(sgn*g)
}

###############################################################

FitMatchDC2 = function(coefs,allcoefs,which,pars,proc,sgn=1)    # Hessian for matching unmeasured components
{
    allcoefs[,which] = matrix(coefs,nrow(allcoefs),length(which))
    
    ind = array(1:length(allcoefs),dim(allcoefs))
    ind = as.vector(ind[,which])

    H = proc$d2fdc2(allcoefs,proc$bvals,pars,proc$more)

    return(as.matrix(sgn*H[ind,ind]))
}

###############################################################

FitMatchList = function(coefs,allcoefs,which,pars,proc,sgn=1)  # Inner Objective
{
  value = FitMatchErr(coefs,allcoefs,which,pars,proc,sgn=1)
  gradient = FitMatchDC(coefs,allcoefs,which,pars,proc,sgn=1)
  hessian = FitMatchDC2(coefs,allcoefs,which,pars,proc,sgn=1)

  return(list(value=value,gradient=gradient,hessian=hessian))

}

###############################################################

ProfileErr.AllPar = function(pars,times,data,coefs,lik,proc,in.meth = "house",control.in=NULL,sgn=1)
{                                                          # Outer Optimization Objective

    coefs = as.vector(coefs)  # First run the inner optimization
    
    if(file.exists('curcoefs.tmp')) {                      # This is enough to make me swear, we'll just
      altcoefs = as.matrix(read.table('curcoefs.tmp'))               # write out to a .tmp file
      if( !(length(altcoefs)==length(coefs)) ){
        stop('Variables in curcoefs.tmp do not conform; file exists from previous experiments?')
       }
    }
    else{
      altcoefs = coefs
    }
    
    if(file.exists('counter.tmp')){                       # And we'll keep track of the
      counter = read.table('counter.tmp')                 # evaluation history
      niter = counter[nrow(counter),1]
    }
    else{
      counter = matrix(c(1,0,pars),1,length(pars)+2)
      niter = 0
    }
    
           
    altdevals = as.matrix(lik$bvals%*%matrix(altcoefs,ncol(lik$bvals),length(altcoefs)/ncol(lik$bvals)))

    colnames(altdevals) = proc$more$names
    
    f1 = SplineCoefsErr(coefs,times,data,lik,proc,pars)
    f2 = SplineCoefsErr(as.vector(altcoefs),times,data,lik,proc,pars)    

    if(f2 < f1){ coefs = altcoefs }

    Ires = inneropt(data,times,pars,coefs,lik,proc,in.meth,control.in)
    ncoefs = Ires$coefs


    devals = as.matrix(lik$bvals%*%ncoefs)
    colnames(devals) = proc$more$names
    
    f = sum(lik$fn(data,times,devals,pars,lik$more))

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

###############################################################

ProfileDP.AllPar = function(pars,times,data,coefs,lik,proc,in.meth=NULL,control.in=NULL,sgn=1,sumlik=1)
{

    if(file.exists('optcoefs.tmp')){
      altcoefs = as.matrix(read.table('optcoefs.tmp'))
      if( !(length(altcoefs)==length(coefs)) ){
        stop('Variables in curcoefs.tmp do not conform; file exists from previous experiments?')
      }else{
        coefs = altcoefs
      }      
    }

    coefs = as.matrix(coefs)
    devals = as.matrix(lik$bvals%*%coefs)
    colnames(devals) = proc$more$names

    coefs = as.vector(coefs)    

    d2Hdc2 = SplineCoefsDC2(coefs,times,data,lik,proc,pars)
    d2Hdcdp = SplineCoefsDCDP(coefs,times,data,lik,proc,pars)
  
    dcdp = -ginv(d2Hdc2)%*%d2Hdcdp
  
   
    if(sumlik){
        dlikdc = as.vector(t(lik$bvals)%*%lik$dfdx(data,times,devals,pars,lik$more))
        df = t(dcdp)%*%(dlikdc) + apply(lik$dfdp(data,times,devals,pars,lik$more),2,sum)    
    
        df = as.vector(sgn*df)
     #   print(df)
        return(df)
    }
    else{
        dlikdx = lik$dfdx(data,times,devals,pars,lik$more)

        dlikdp = lik$dfdp(data,times,devals,pars,lik$more)

        dlikdc = c()
        for(i in 1:dim(dlikdx)[2]){
            dlikdc = cbind(dlikdc,as.matrix(diag(dlikdx[,i])%*%lik$bvals))            
        }       
        
        df = dlikdc%*%dcdp + dlikdp
        
        return(df)
    }
}

###############################################################

ProfileErr = function(pars,allpars,times,data,coefs,lik,proc,in.meth = "house",
        control.in=NULL,sgn=1,active=1:length(allpars))
{
    allpars[active] = pars
    
    f = ProfileErr.AllPar(allpars,times,data,coefs,lik,proc,in.meth,control.in,sgn)
    return(f)
}

###############################################################

ProfileDP = function(pars,allpars,times,data,coefs,lik,proc,in.meth = "house",
                control.in=NULL,sgn=1,sumlik=1,active=1:length(allpars))
{
    allpars[active] = pars

    g = ProfileDP.AllPar(allpars,times,data,coefs,lik,proc,in.meth,control.in,sgn,sumlik)


    if(sumlik){ 
        names(g) = names(allpars)
        g = g[active] 
    }
    else{ 
        colnames(g) = names(allpars)
        g = g[,active] 
    }
    
    return(g)        
}

###############################################################

ProfileList = function(pars,allpars,times,data,coefs,lik,proc,in.meth = "house",
        control.in=NULL,sgn=1,active=1:length(allpars))
{
  value = ProfileErr(pars,allpars,times,data,coefs,lik,proc,in.meth,control.in,sgn,active)
  gradient = ProfileDP(pars,allpars,times,data,coefs,lik,proc,in.meth,control.in,sgn,active)
  
  return(list(value=value,gradient=gradient))
}

###############################################################

ProfileSSE = function(pars,allpars,times,data,coefs,lik,proc,in.meth='nlminb',control.in=NULL,   # Squared Error
                    active=1:length(pars),dcdp=NULL,oldpars=NULL,use.nls=TRUE,sgn=1)                            # outer criterion
{
  allpars[active] = pars

  f = ProfileSSE.AllPar(pars=allpars,times=times,data=data,coefs=coefs,lik=lik,proc=proc,in.meth=in.meth,
      control.in=control.in,dcdp=dcdp,oldpars=oldpars,use.nls=use.nls,sgn=sgn)
  
    if(use.nls){
#        attr(f,'gradient') = attr(f,'gradient')[,active]
    }
    else{
        f$gradient = f$gradient[,active]
        f$dcdp = f$dcdp[,active]
    }
    return(f)
}

###############################################################

ProfileSSE.AllPar = function(pars,times,data,coefs,lik,proc,in.meth='nlminb',control.in=NULL,   # Squared Error
                    dcdp=NULL,oldpars=NULL,use.nls=TRUE,sgn=1)                            # outer criterion
{           

    coefs = as.vector(coefs)  # First run the inner optimization
    f1 = SplineCoefsErr(coefs,times,data,lik,proc,pars)

    if(use.nls){                                        # If using NLS, we need to keep track
        if(file.exists('curcoefs.tmp')){
          altcoefs = as.matrix(read.table('curcoefs.tmp'))
          if( !(length(altcoefs)==length(coefs)) ){
            stop('Variables in curcoefs.tmp do not conform; file exists from previous experiments?')
          }
        }
        else{
          altcoefs = coefs
        }
        if(file.exists('counter.tmp')){                       # And we'll keep track of the
          counter = read.table('counter.tmp')                 # evaluation history
          niter = counter[nrow(counter),1]
        }
        else{
          counter = matrix(c(1,0,pars),1,length(pars)+2)
          niter = 0
        }

        f2 = SplineCoefsErr(as.vector(altcoefs),times,data,lik,proc,pars)    
    
        if(f2 < f1){
          coefs = altcoefs
          f1 = f2
        }
        
        altdevals = as.matrix(lik$bvals%*%matrix(altcoefs,ncol(lik$bvals),length(altcoefs)/ncol(lik$bvals)))
        colnames(altdevals) = proc$more$names
    }

    if( !is.null(dcdp) ){
        tcoefs = coefs + dcdp%*%(pars-oldpars);
        f2 = SplineCoefsErr(as.vector(tcoefs),times,data,lik,proc,pars)
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
    f = as.vector( as.matrix(data - lik$more$fn(times,devals,pars,lik$more$more))*sqrt(weights) )
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

    d2Hdc2 = SplineCoefsDC2(ncoefs,times,data,lik,proc,pars)
    d2Hdcdp = SplineCoefsDCDP(ncoefs,times,data,lik,proc,pars)

    dcdp = ginv(d2Hdc2)%*%d2Hdcdp

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
            write.table(ncoefs,file='curcoefs.tmp',row.names=FALSE,col.names=FALSE)
            niter = counter[nrow(counter),1]
        }
        
        if(niter==0){ 
          counter[1,2] = tf
          write.table(counter,file='counter.tmp',col.names=FALSE,row.names=FALSE) 
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

###############################################################

Profile.GausNewt = function(pars,times,data,coefs,lik,proc,   # Gauss-Newton routine for squared error outer
    in.meth='nlminb',control.in=NULL,active=1:length(pars),   # objective
    control=list(reltol=1e-6,maxit=50,maxtry=15,trace=1)) 
{
    res0 = ProfileSSE.AllPar(pars,times,data,coefs,lik,proc,in.meth,control.in,use.nls=FALSE)

    ind = names(pars)
    if(is.null(ind)){ ind = 1:length(pars) }                # Allows some pars to be
    ind = ind[active]                                       # fixed

    pars0 = pars
    pars1 = pars
    F0 = sum(res0$value^2)
    F1 = F0
    iter = 0
    fundif = 1

    gradnorm0 = mean(abs(t(res0$gradient[,ind])%*%res0$value))
    gradnorm1 = 1
    dcdp = c()

    while( gradnorm1 > control$reltol & fundif > control$reltol & iter < control$maxit){
  
        iter = iter + 1

        Dpars = solve(t(res0$gradient[,ind])%*%res0$gradient[,ind],t(res0$gradient[,ind])%*%res0$value)
    
        ntry = 0
        
        while(F1 >= F0 & t(Dpars)%*%Dpars > control$reltol & ntry < control$maxtry){
            pars1[ind] = pars0[ind] - 0.5*Dpars
            tres = ProfileSSE.AllPar(pars1,times,data,res0$coefs,lik,proc,in.meth,control.in,res0$dcdp,pars0,use.nls=FALSE)
            F1 = sum(tres$value^2)
#            print(c(F0,F1,t(Dpars)%*%Dpars,pars1))
            Dpars = Dpars/2
            ntry = ntry + 1
        }

        res0 = tres
        dcdp = tres$dcdp
        gradnorm1 = abs(mean(t(res0$gradient[,ind])%*%res0$value))
        fundif = (F0-F1)/abs(F0)

        pars0 = pars1
        res0 = tres
        F0 = F1
        gradnorm0 = gradnorm1
        
        if(control$trace>0){ print(c(iter,ntry,F0,pars0)) }
    }

    newpars = pars0
     
    return(list(pars=newpars,in.res=res0,value=F0))
}

###############################################################

SplineEst.NewtRaph = function(coefs,times,data,lik,proc,pars,   # Newton-Raphson routine for 
    control=list(reltol=1e-12,maxit=1000,maxtry=10,trace=0))    # Inner optimization
{
    f0 = SplineCoefsErr(coefs,times,data,lik,proc,pars)
    g = SplineCoefsDC(coefs,times,data,lik,proc,pars)
    H = SplineCoefsDC2(coefs,times,data,lik,proc,pars)

    DC = -ginv(H)%*%g

    gradnorm1 = 1
    fundif = 1
    iter = 0
    f1 = f0

    coefs0 = as.vector(coefs)

    while( gradnorm1 > control$reltol & fundif > 0 & iter < control$maxit){
  
        iter = iter + 1    
        DC = -ginv(H)%*%g
   
        ntry = 0
        
        coefs1 = coefs0 

        while( (f1 >= f0) & (t(DC)%*%DC > control$reltol) & (ntry < control$maxtry) ){
            
            
            coefs1 = coefs0 + DC
            f1 = SplineCoefsErr(coefs1,times,data,lik,proc,pars)
            DC = DC/2
            ntry = ntry + 1
  #          print(c(f1,f0, t(DC)%*%DC,control$reltonl,ntry,control$maxtry))
  #          print((f1 >= f0) & (t(DC)%*%DC > control$reltol) & (ntry < control$maxtry))
            
        }

        coefs0 = coefs1
        
        g = SplineCoefsDC(coefs0,times,data,lik,proc,pars)
        H = SplineCoefsDC2(coefs0,times,data,lik,proc,pars)    

        gradnorm1 = mean(abs(DC))
        fundif = (f0-f1)/abs(f0)
  
        f0 = f1      
        if(control$trace>1){ print(c(iter,ntry,f0,gradnorm1,fundif)) }

    }
    if(control$trace>0){ print(c(f0,gradnorm1,iter)) }
     
    return(list(coefs=coefs0,g=g,value=f0,H=H))
}

###############################################################

NeweyWest.Var = function(H,g,maxlag)       
{
    V = solve(H)

    I = 0*H
    
    if(is.null(maxlag)){ 
        n = nrow(g)
        maxlag = max(5,n^(0.25))
    }

    if(maxlag > 0){
        for(i in 1:(ncol(g)-1)){
            for(j in (i+1):ncol(g)){
                I[i,j] = Newey.West(g[,i],g[,j],maxlag)
                I[j,i] = I[i,j]
            }    
        }
     
    }

    return( V%*%(I+ t(g)%*%g)%*%V )

}

###############################################################

Profile.covariance <- function(pars,active=NULL,times,data,coefs,lik,proc,in.meth='nlminb',control.in=NULL,eps=1e-6)
{
    if(is.null(active)){ active = 1:length(pars) }

    apars = pars[active]
    
    H = matrix(0,length(apars),length(apars))

    g = ProfileDP(pars=apars,allpars=pars,times=times,data=data,coefs=coefs,lik=lik,proc=proc,active=active,sumlik=FALSE)
    gg = apply(g,2,sum)

    for(i in 1:length(apars)){        
      if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
      if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
      if(file.exists('counter.tmp')){file.remove('counter.tmp')}

        tpars = apars
        tpars[i] = tpars[i] + eps

        tf = ProfileErr(tpars,pars,times,data,coefs,lik,proc,in.meth=in.meth,control.in=control.in,active=active)  
        tg = ProfileDP(tpars,pars,times,data,coefs,lik,proc,active=active,sumlik=TRUE)                            

        H[,i] = (tg-gg)/eps

      if(file.exists('curcoefs.tmp')){file.remove('curcoefs.tmp')}
      if(file.exists('optcoefs.tmp')){file.remove('optcoefs.tmp')}
      if(file.exists('counter.tmp')){file.remove('counter.tmp')}
    }
               
    Covar = NeweyWest.Var( 0.5*(t(H)+H) ,g,5)
        
    return( Covar )
}

############################################################################################
#
# Some Utilities
#
############################################################################################


#blocks2mat = function(H)  # List of matrices -> large matrix
#{
#
#    out = c()
#
#    for(i in 1:dim(H)[3]){
#
#    tout = c()
#    for(j in 1:dim(H)[4]){
#        tout = cbind(tout,H[,,i,j])
#    }
#
#    out = rbind(out,tout)
#    }
#
#}

blocks2mat = function(H)  # List of matrices -> large matrix
{

    out = c()

    for(i in 1:length(H)){
      tout = H[[i]][[1]]
      if(length(H[[i]])>1){
        for(j in 2:length(H[[i]])){
            if(inherits(H[[i]][[j]],'dgCMatrix')|inherits(H[[i]][[j]],'dgeMatrix')){ tout=cBind(tout,H[[i]][[j]]) }
            else{ tout = cbind(tout,H[[i]][[j]]) }
        }
      }
      if(i > 1){ 
        if(inherits(tout,'dgCMatrix')|inherits(tout,'dgeMatrix')){ out=rBind(out,tout) }
        else{ out = rbind(out,tout) }
      } 
      else{ out = tout }
    }

    return(out)
}

## Newey West Calculations, with thanks to Steve Ellner


## GAUSS trimr function: trims n1 rows from the start and n2 rows from the end
## of a matrix or vector 
trimr <- function (a,n1,n2) {
        da<-dim(a); 
        if(is.null(da)) {a[(n1+1):(length(a)-n2)]}
        else {a[(n1+1):(da[1]-n2),]}
}

Newey.West <-function(x,y,maxlag) {
        w=1-(1:maxlag)/(maxlag+1); w=w/length(x); 
        out=mean(x*y); 
        for(i in 1:maxlag) {
            out=out+w[i]*sum(trimr(x,i,0)*trimr(y,0,i))+w[i]*sum(trimr(y,i,0)*trimr(x,0,i))
        }
        return(out)     
} 


