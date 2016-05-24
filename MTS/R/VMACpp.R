#### VMA programs
"VMACpp" <- function(da,q=1,include.mean=T,fixed=NULL,beta=NULL,sebeta=NULL,prelim=F,details=F,thres=2.0){
   # Estimation of a vector MA model using conditional MLE (Gaussian dist)
   #
   # April 18: add subcommand "prelim" to see simplification after the AR approximation.
   # When prelim=TRUE, fixed is assigned based on the results of AR approximation.
   # Here "thres" is used only when prelim = TRUE.
   ##
   ### Create the "mFilter" program to simplify computation of residuals. April 8, 2012.
   #
   if(!is.matrix(da))da=as.matrix(da)
   nT=dim(da)[1]
   k=dim(da)[2]
   if(q < 1)q=1
   kq=k*q
#
 THini <- function(y,x,q,include.mean){
   # use residuals of a long VAR model to obtain initial estimates of
   # VMA coefficients.
   if(!is.matrix(y))y=as.matrix(y)
   if(!is.matrix(x))x=as.matrix(x)
   nT=dim(y)[1]
   k=dim(y)[2]
   ist=1+q
   ne=nT-q
   if(include.mean){
      xmtx=matrix(1,ne,1)
   }
   else {
      xmtx=NULL
   }
   ymtx=y[ist:nT,]
   for (j in 1:q){
      xmtx=cbind(xmtx,x[(ist-j):(nT-j),])
   }
   xtx=crossprod(xmtx,xmtx)
   xty=crossprod(xmtx,ymtx)
   xtxinv=solve(xtx)
   beta=xtxinv%*%xty
   resi= ymtx - xmtx%*%beta
   sse=crossprod(resi,resi)/ne
   dd=diag(xtxinv)
   sebeta=NULL
   for (j in 1:k){
      se=sqrt(dd*sse[j,j])
      sebeta=cbind(sebeta,se)
   }
    THini <- list(estimates=beta,se=sebeta)
  }

   if(length(fixed) < 1){ 
    m1=VARorder(da,q+12,output=FALSE)
    porder=m1$aicor
    if(porder < 1)porder=1
    m2=VAR(da,porder,output=FALSE)
    y=da[(porder+1):nT,]
    x=m2$residuals
    m3=THini(y,x,q,include.mean)
    beta=m3$estimates
    sebeta=m3$se
    nr=dim(beta)[1]
   ### Preliminary simplification
   if(prelim){
      fixed = matrix(0,nr,k)
      for (j in 1:k){
         tt=beta[,j]/sebeta[,j]
         idx=c(1:nr)[abs(tt) >= thres]
         fixed[idx,j]=1
       }
     }
   #
    if(length(fixed) < 1){fixed=matrix(1,nr,k)}
   }
   else{
    nr=dim(beta)[1]
   }
   #
   par=NULL
   separ=NULL
   fix1=fixed
   #
   #
   VMAcnt=0
   ist=0
   if(include.mean){
      jdx=c(1:k)[fix1[1,]==1]
      VMAcnt=length(jdx)
      if(VMAcnt > 0){
         par=beta[1,jdx]
         separ=sebeta[1,jdx]
      }
      TH=-beta[2:(kq+1),]
      seTH=sebeta[2:(kq+1),]
      ist=1
   }
   else {
      TH=-beta
      seTH=sebeta
   }
   #########
   for (j in 1:k){
      idx=c(1:(nr-ist))[fix1[(ist+1):nr,j]==1]
      if(length(idx) > 0){
         par=c(par,TH[idx,j])
         separ=c(separ,seTH[idx,j])
      }
   }
   #
   ParMA <- par
  
  LLKVMACpp <- function(par,zt=zt,q=q,fixed=fix1,include.mean=include.mean){
   k=ncol(zt)
   nT=nrow(zt)
   kq=k*q
   fix <- fixed
   
   ListObs = .Call("GetVMAObs", zt, fix1, par, q, include.mean)
   zt  = do.call(rbind, ListObs)
   
   #### Slow!
   ListTH = .Call("GetVMATH", zt, fix1, par, q, include.mean)
   TH  = do.call(rbind, ListTH)
       
   mm=eigen(TH)
   V1=mm$values
   P1=mm$vectors
   v1=Mod(V1)
   ich=0
   for (i in 1:kq){
      if(v1[i] > 1)V1[i]=1/V1[i]
      ich=1
   }
   if(ich > 0){
      ###cat("Invertibility checked and adjusted: ","\n")
      P1i=solve(P1)
      GG=diag(V1)
      TH=Re(P1%*%GG%*%P1i)
      Theta=t(TH[1:k,])
      ##cat("adjusted Theta","\n")
      ##print(TH[1:k,])
      ### re-adjust the MA parameter
      icnt <- VMAcnt
      ist=0
      if(icnt > 0)ist=1
      for (j in 1:k){
         idx=c(1:kq)[fix[(ist+1):(ist+kq),j]==1]
         jcnt=length(idx)
         if(jcnt > 0){
            par[(icnt+1):(icnt+jcnt)]=TH[j,idx]
            icnt=icnt+jcnt
         }
      }
      ##
      ParMA <- par
   }
   
   ##
   at=mFilter(zt,t(Theta))
   #
   sig=t(at)%*%at/nT
   ##ll=dmnorm(at,rep(0,k),sig)
   ll=dmvnorm(at,rep(0,k),sig)
   LLKVMACpp =-sum(log(ll))
   LLKVMACpp
}



  
   #
   cat("Number of parameters: ",length(par),"\n")
   cat("initial estimates: ",round(par,4),"\n")
   ### Set up lower and upper bounds
   lowerBounds=par; upperBounds=par
   npar=length(par)
   mult=2.0
   if((npar > 10)||(q > 2))mult=1.2
   for (j in 1:npar){
      lowerBounds[j] = par[j]-mult*separ[j]
      upperBounds[j] = par[j]+mult*separ[j]
   }
   cat("Par. Lower-bounds: ",round(lowerBounds,4),"\n")
   cat("Par. Upper-bounds: ",round(upperBounds,4),"\n")
   ###mm=optim(par,LLKvma,method=c("L-BFGS-B"),lower=lowerBounds,upper=upperBounds,hessian=TRUE)
   ###mm=optim(par,LLKvma,method=c("BFGS"),hessian=TRUE)
   ##est=mm$par
   ##H=mm$hessian
   # Step 5: Estimate Parameters and Compute Numerically Hessian:
   if(details){
      fit = nlminb(start = ParMA, objective = LLKVMACpp,zt=da,fixed=fixed,include.mean=include.mean,q=q,
      lower = lowerBounds, upper = upperBounds, control = list(trace=3))
   }
   else {
      fit = nlminb(start = ParMA, objective = LLKVMACpp, zt=da, fixed=fixed, include.mean=include.mean, q=q, 
       lower = lowerBounds, upper = upperBounds)
   }
   epsilon = 0.0001 * fit$par
   npar=length(par)
   Hessian = matrix(0, ncol = npar, nrow = npar)
   for (i in 1:npar) {
      for (j in 1:npar) {
         x1 = x2 = x3 = x4  = fit$par
         x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
         x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
         x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
         x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
         Hessian[i, j] = (LLKVMACpp(x1,zt=da,q=q,fixed=fixed,include.mean=include.mean)
                          -LLKVMACpp(x2,zt=da,q=q,fixed=fixed,include.mean=include.mean)
                          -LLKVMACpp(x3,zt=da,q=q,fixed=fixed,include.mean=include.mean)
                          +LLKVMACpp(x4,zt=da,q=q,fixed=fixed,include.mean=include.mean))/
         (4*epsilon[i]*epsilon[j])
      }
   }
   est=fit$par
   cat("Final   Estimates: ",est,"\n")
   # Step 6: Create and Print Summary Report:
   se.coef = sqrt(diag(solve(Hessian)))
   tval = fit$par/se.coef
   matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
   dimnames(matcoef) = list(names(tval), c(" Estimate",
   " Std. Error", " t value", "Pr(>|t|)"))
   cat("\nCoefficient(s):\n")
   printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
   #
   ### recover to the format of unconstrained case for printing purpose.
   cat("---","\n")
   cat("Estimates in matrix form:","\n")
   icnt=0
   ist=0
   cnt=NULL
   if(include.mean){
      ist=1
      cnt=rep(0,k)
      secnt=rep(1,k)
      jdx=c(1:k)[fix1[1,]==1]
      icnt=length(jdx)
      if(icnt > 0){
         cnt[jdx]=est[1:icnt]
         secnt[jdx]=se.coef[1:icnt]
         cat("Constant term: ","\n")
         cat("Estimates: ",cnt,"\n")
      }
   }
   cat("MA coefficient matrix","\n")
   TH=matrix(0,kq,k)
   seTH=matrix(1,kq,k)
   for (j in 1:k){
      idx=c(1:kq)[fix1[(ist+1):nr,j]==1]
      jcnt=length(idx)
      if(jcnt > 0){
         TH[idx,j]=est[(icnt+1):(icnt+jcnt)]
         seTH[idx,j]=se.coef[(icnt+1):(icnt+jcnt)]
         icnt=icnt+jcnt
      }
   }
   icnt=0
   for (i in 1:q){
      cat("MA(",i,")-matrix","\n")
      theta=t(TH[(icnt+1):(icnt+k),])
      print(theta,digits=3)
      icnt=icnt+k
   }
   ## Compute the residuals
   zt=da
   if(include.mean){
      for (i in 1:k){
         zt[,i]=zt[,i]-cnt[i]
      }
   }
   ### Use mFilter to compute residuals (April 18, 2012)
   at=mFilter(zt,t(TH))
   sig=t(at)%*%at/nT
   cat(" ","\n")
   cat("Residuals cov-matrix:","\n")
   print(sig)
   dd=det(sig)
   d1=log(dd)
   aic=d1+2*npar/nT
   bic=d1+log(nT)*npar/nT
   cat("----","\n")
   cat("aic= ",aic,"\n")
   cat("bic= ",bic,"\n")
   ### prepare for output storage
   Theta=t(TH)
   if(include.mean){
      TH=rbind(cnt,TH)
      seTH=rbind(secnt,seTH)
   }
   
   VMACpp <- list(data=da,MAorder=q,cnst=include.mean,coef=TH,secoef=seTH,residuals=at,Sigma=sig,Theta=Theta,mu=cnt,aic=aic,bic=bic)
}