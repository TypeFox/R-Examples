### An MTS package by Ruey S. Tsay
##library(mvtnorm)
###
"MTSplot" <- function(data,caltime=NULL){
   ## plot the multivariate time series
   ### caltime: calendar time
   if(!is.matrix(data))data=as.matrix(data)
   if(is.ts(data)){
      plot(data)
   }
   else{
      nT=dim(data)[1]
      tdx=c(1:nT)
      if(length(caltime) > 1)tdx=caltime
      k=dim(data)[2]
      if(k < 4){
       par(mfcol=c(k,1))
       for (j in 1:k){
        plot(tdx,data[,j],xlab='time',ylab=colnames(data)[j],type='l')
        }
       }
      if(k == 4){
         par(mfcol=c(2,2))
         for (j in 1:k){
          plot(tdx,data[,j],xlab='time',ylab=colnames(data)[j],type='l')
          }
      }
      if((k > 4) && (k < 13)){
         par(mfcol=c(3,2),mai=c(0.3,0.3,0.3,0.3))
         k1=6
        jcnt=0
        for (j in 1:k){
         plot(tdx,data[,j],xlab='time',ylab=colnames(data)[j],type='l',cex.axis=0.8)
         jcnt=jcnt+1
         if((jcnt == k1) && (k > 6)){
            jcnt=0
            cat("Hit return for more plots: ","\n")
            readline()
         }
      }
    }
    if(k > 12){
     par(mfcol=c(1,1))
     yl=range(data)*1.05
     plot(tdx,data[,1],xlab='time',ylab=' ',type='l',ylim=yl)
     for (j in 2:k){
      lines(tdx,data[,j],lty=j,col=j)
      }
    }
   #end of the program
   }
   par(mfcol=c(1,1))
}

####
"VAR" <- function(x,p=1,output=T,include.mean=T,fixed=NULL){
   # Fits a vector AR(p) model, computes AIC, and residuals
   # fixed[i,j] = 1 denotes the parameter needs estimation, = 0, means fixed to 0.
   if(!is.matrix(x))x=as.matrix(x)
   Tn=dim(x)[1]
   k=dim(x)[2]
   if(p < 1)p=1
   idm=k*p
   ne=Tn-p
   ist=p+1
   y=x[ist:Tn,]
   if(include.mean){
      idm=idm+1
      xmtx=cbind(rep(1,ne),x[p:(Tn-1),])
   }
   else {
      xmtx=x[p:(Tn-1),]
   }
   if(p > 1){
      for (i in 2:p){
         xmtx=cbind(xmtx,x[(ist-i):(Tn-i),])
      }
   }
   #
   ndim=ncol(xmtx)
   if(length(fixed)==0){
      paridx=matrix(1,ndim,k)
   }
   else {
      paridx=fixed
   }
   #perform estimation component-by-component
   res=NULL
   beta=matrix(0,ndim,k)
   sdbeta=matrix(0,ndim,k)
   npar=0
   for (i in 1:k){
      idx=c(1:ndim)[paridx[,i]==1]
      resi=y[,i]
      if(length(idx)> 0){
         xm=as.matrix(xmtx[,idx])
         npar=npar+dim(xm)[2]
         xpx=t(xm)%*%xm
         xpxinv=solve(xpx)
         xpy=t(xm)%*%as.matrix(y[,i],ne,1)
         betai=xpxinv%*%xpy
         beta[idx,i]=betai
         resi=y[,i]-xm%*%betai
         nee=dim(xm)[2]
         sse=sum(resi*resi)/(Tn-p-nee)
         dd=diag(xpxinv)
         sdbeta[idx,i]=sqrt(dd*sse)
      }
      res=cbind(res,resi)
   }
   sse=t(res)%*%res/(Tn-p)
   ###cat("npar =",npar,"\n")
   #
   aic=0
   bic=0
   hq=0
   Phi=NULL
   Ph0=NULL
   #
   jst=0
   if(include.mean) {
      Ph0=beta[1,]
      se=sdbeta[1,]
      if(output){
         cat("Constant term:","\n")
         cat("Estimates: ",Ph0,"\n")
         cat("Std.Error: ",se,"\n")
      }
      jst=1
   }
   ### adjustment npar for computing information criterion
   if(include.mean){
      for (i in 1:k){
         if(abs(Ph0[i]) > 0.00000001)npar=npar-1
      }
   }
   ####cat("adjusted npar = ",npar,"\n")
   if(output)cat("AR coefficient matrix","\n")
   ### Phi is a storage for AR coefficient matrices
   for (i in 1:p){
      phi=t(beta[(jst+1):(jst+k),])
      se=t(sdbeta[(jst+1):(jst+k),])
      if(output){
         cat("AR(",i,")-matrix","\n")
         print(phi,digits=3)
         cat("standard error","\n")
         print(se,digits=3)
      }
      jst=jst+k
      Phi=cbind(Phi,phi)
   }
   if(output){
      cat(" ","\n")
      cat("Residuals cov-mtx:","\n")
      print(sse)
      cat(" ","\n")
   }
   dd=det(sse)
   d1=log(dd)
   aic=d1+(2*npar)/Tn
   bic=d1+log(Tn)*npar/Tn
   hq=d1+2*log(log(Tn))*npar/Tn
   if(output){
      cat("det(SSE) = ",dd,"\n")
      cat("AIC = ",aic,"\n")
      cat("BIC = ",bic,"\n")
      cat("HQ  = ",hq,"\n")
      # end of if(output)
   }
   
   VAR<-list(data=x,cnst=include.mean,order=p,coef=beta,aic=aic,bic=bic,hq=hq,residuals=res,secoef=sdbeta,Sigma=sse,Phi=Phi,Ph0=Ph0)
}

#####
"refVAR" <- function(model,fixed=NULL,thres=1.0){
   # This program automatically refines the "model" by removing estimates
   # with abs(t-ration) < thres.
   #
   # model is an VAR output object.
   #
   x=as.matrix(model$data)
   nT=dim(x)[1]
   k=dim(x)[2]
   p=model$order
   if(p < 1)p=1
   cnst=model$cnst
   fix=fixed
   if(length(fixed)== 0){
      coef=as.matrix(model$coef)
      secoef=as.matrix(model$secoef)
      nr=dim(coef)[1]
      nc=dim(coef)[2]
      for (i in 1:nr){
         for (j in 1:nc){
            if(secoef[i,j] < 10^(-8))secoef[i,j]=1
         }
      }
      fix=matrix(1,nr,k)
      # use Backward elimination to simplify the model: equation by equation
      ### First: setup the regressor matrix
      xmtx=NULL
      ist=p+1
      y=x[ist:nT,]
      ne=nT-p
      if(cnst)xmtx=matrix(1,ne,1)
      for (j in 1:p){
         xmtx=cbind(xmtx,x[(ist-j):(nT-j),])
      }
      xmtx=as.matrix(xmtx)
      ### perform first elimination based on the previous estimation
      for(j in 1:k){
         tt=abs(coef[,j]/secoef[,j])
         idx=c(1:nr)[tt == min(tt)]
         idx1=idx[1]
         if(tt[idx1] < thres)fix[idx,j]=0
      }
      ### Perform further elimination
      for (j in 1:k){
         npar=sum(fix[,j])
         while(npar > 0){
            jdx=c(1:nr)[fix[,j]==1]
            xp=as.matrix(xmtx[,jdx])
            nxp=dim(xp)[2]
            m1=lm(y[,j]~-1+xp)
            m2=summary(m1)
            est=m1$coefficients
            se1=sqrt(diag(m2$cov.unscaled))*m2$sigma
            tt=abs(est/se1)
            idx=c(1:nxp)[tt == min(tt)]
            idx1=idx[1]
            if(tt[idx1] < thres){
               fix[jdx[idx],j]=0
               npar=sum(fix[,j])
            }
            else {
               npar=0
            }
            ### end of while-statement
         }
         ## end of the j-loop
      }
      ## end of the if(length(fixed)==0) statement
   }
   
   mm=VAR(x,p,output=T,include.mean=cnst,fixed=fix)
   
   refVAR <- list(data=mm$data,order=p,cnst=cnst,coef=mm$coef,aic=mm$aic,bic=mm$bic,hq=mm$hq,residuals=mm$residuals,secoef=mm$secoef,Sigma=mm$Sigma,Phi=mm$Phi,Ph0=mm$Ph0)
}

######
"VARs" <- function(x,lags,include.mean=T,output=T,fixed=NULL){
   # Fits a vector AR model with selected lags, computes its AIC, BIC, and residuals
   # Created on March 30, 2009 by Ruey S. Tsay
   #
   if(!is.matrix(x))x=as.matrix(x)
   nT=dim(x)[1]; k=dim(x)[2]
   nlags=length(lags)
   #
   if(nlags < 0){
      lags=c(1)
      nlags=1
   }
   #
   lags=sort(lags)
   idm=k*nlags+1
   p=lags[nlags]
   if(p < 1)p=1
   ne=nT-p
   ist=p+1
   y=x[ist:nT,]
   jj=lags[1]
   if(include.mean){
      xmtx=cbind(rep(1,ne),x[(ist-jj):(nT-jj),])
   }
   else {
      xmtx=x[(ist-jj):(nT-jj),]
   }
   if(nlags > 1){
      for (i in 2:nlags){
         jj=lags[i]
         xmtx=cbind(xmtx,x[(ist-jj):(nT-jj),])
      }
   }
   xmtx=as.matrix(xmtx)
   ndim=dim(xmtx)[2]
   #
   if(length(fixed)==0){
      paridx=matrix(1,ndim,k)
   }
   else {
      paridx=fixed
   }
   #perform estimation component-by-component
   res=NULL
   beta=matrix(0,ndim,k)
   sdbeta=matrix(0,ndim,k)
   for (i in 1:k){
      idx=c(1:ndim)[paridx[,i]==1]
      resi=y[,i]
      if(length(idx) > 0){
         xm=as.matrix(xmtx[,idx])
         xpx = t(xm)%*%xm
         xpxinv=solve(xpx)
         xpy=t(xm)%*%as.matrix(y[,i],ne,1)
         betai=xpxinv%*%xpy
         beta[idx,i]=betai
         resi=y[,i]-xm%*%betai
         sse=sum(resi*resi)/nT
         dd=diag(xpxinv)
         sdbeta[idx,i]=sqrt(sse*dd)
      }
      res=cbind(res,resi)
   }
   sse=t(res)%*%res/nT
   Phi=NULL
   Ph0=rep(0,k)
   if(output){
      jst=0
      if(include.mean){
         Ph0=beta[1,]
         se=sdbeta[1,]
         cat("Constant term:","\n")
         print(Ph0,digits=4)
         cat("Std error:","\n")
         print(se,digits=4)
         jst=1
      }
      cat("AR coefficient matrix:","\n")
      ### Phi is a storage for the AR coefficient matrices
      Phi=matrix(0,k,p*k)
      for (i in 1:nlags){
         ord=lags[i]
         cat("AR(",ord,")-matrix:","\n")
         phi=t(beta[(jst+1):(jst+k),])
         se=t(sdbeta[(jst+1):(jst+k),])
         print(phi,digits=3)
         cat("Standard error:","\n")
         print(se,digits=3)
         jst=jst+k
         cat("      ","\n")
         kdx=(ord-1)*k
         Phi[,(kdx+1):(kdx+k)]=phi
      }
      cat("Residuals cov-mtx:","\n")
      print(sse)
      
      cat("       ","\n")
      dd=det(sse)
      cat("det(SSE) = ",dd,"\n")
      d1=log(dd)
      aic=d1+(2*nlags*k*k)/nT
      bic=d1+log(nT)*nlags*k*k/nT
      cat("AIC = ",aic,"\n")
      cat("BIC = ",bic,"\n")
      # end of if(output)
   }
   
   VARs<-list(data=x,lags=lags,order=p,cnst=include.mean,coef=beta,aic=aic,bic=bic,residuals=res,secoef=sdbeta,Sigma=sse,Phi=Phi,Ph0=Ph0)
   
}

#####
"refVARs" <- function(model,fixed=NULL,thres=1.0){
   #This program either (1) automatically refines a fitted VARs model by setting
   # parameters with abs(t-ratio) < thres to zero, or uses manually specified "fixed" to
   # simplied the VARs model.
   #
   # model: an VARs output object.
   #
   x=as.matrix(model$data)
   nT=dim(x)[1]; k=dim(x)[2]; p=model$order
   lags=sort(model$lags)
   nlags = length(lags)
   cnst=model$cnst
   fix=fixed
   if(length(fixed)==0){
      p=lags[nlags]
      coef=as.matrix(model$coef)
      secoef=as.matrix(model$secoef)
      nr=dim(coef)[1]
      fix=matrix(1,nr,k)
      nc=dim(coef)[2]
      for (i in 1:nr){
         for (j in 1:nc){
            if(secoef[i,j] < 10^(-8))secoef[i,j]=1.0
         }
      }
      # use Backward elimination to simplify the model: equation by equation
      ### First: setup the regressor matrix
      xmtx=NULL
      ist=p+1
      y=x[ist:nT,]
      ne=nT-p
      if(cnst)xmtx=cbind(xmtx,rep(1,ne))
      for (j in 1:nlags){
         xmtx=cbind(xmtx,x[(ist-lags[j]):(nT-lags[j]),])
      }
      xmtx=as.matrix(xmtx)
      ### perform first elimination based on the previous estimation
      for(j in 1:k){
         tt=abs(coef[,j]/secoef[,j])
         idx=c(1:nr)[tt == min(tt)]
         if(tt[idx] < thres)fix[idx,j]=0
      }
      ### Perform further elimination
      for (j in 1:k){
         npar=sum(fix[,j])
         while(npar > 0){
            jdx=c(1:nr)[fix[,j]==1]
            xp=as.matrix(xmtx[,jdx])
            nxp=dim(xp)[2]
            m1=lm(y[,j]~-1+xp)
            m2=summary(m1)
            est=m1$coefficients
            se1=sqrt(diag(m2$cov.unscaled))*m2$sigma
            tt=abs(est/se1)
            idx=c(1:nxp)[tt == min(tt)]
            if(tt[idx] < thres){
               fix[jdx[idx],j]=0
               npar=sum(fix[,j])
            }
            else {
               npar=0
            }
            ### end of while-statement
         }
         ## end of the j-loop
      }
      ## end of the if(length(fixed)==0) statement
   }
   
   mm=VARs(x,lags,include.mean=cnst,output=T,fixed=fix)
   
   refVARs <- list(data=x,lags=lags,cnst=cnst,coef=mm$coef,aic=mm$aic,bic=mm$bic,residuals=mm$residuals,secoef=mm$secoef,Sigma=mm$Sigma,Phi=mm$Phi,Ph0=mm$Ph0)
}


#####
"ccm" <- function(x,lags=12,level=FALSE,output=T){
   # Compute and plot the cross-correlation matrices.
   # lags: number of lags used.
   # level: logical unit for printing.
   #
   if(!is.matrix(x))x=as.matrix(x)
   nT=dim(x)[1]; k=dim(x)[2]
   if(lags < 1)lags=1
   # remove the sample means
   y=scale(x,center=TRUE,scale=FALSE)
   V1=cov(y)
   if(output){
      print("Covariance matrix:")
      print(V1,digits=3)
   }
   se=sqrt(diag(V1))
   SD=diag(1/se)
   S0=SD%*%V1%*%SD
   ## S0 used later
   ksq=k*k
   wk=matrix(0,ksq,(lags+1))
   wk[,1]=c(S0)
   j=0
   if(output){
      cat("CCM at lag: ",j,"\n")
      print(S0,digits=3)
      cat("Simplified matrix:","\n")
   }
   y=y%*%SD
   crit=2.0/sqrt(nT)
   for (j in 1:lags){
      y1=y[1:(nT-j),]
      y2=y[(j+1):nT,]
      Sj=t(y2)%*%y1/nT
      Smtx=matrix(".",k,k)
      for (ii in 1:k){
         for (jj in 1:k){
            if(Sj[ii,jj] > crit)Smtx[ii,jj]="+"
            if(Sj[ii,jj] < -crit)Smtx[ii,jj]="-"
         }
      }
      #
      if(output){
         cat("CCM at lag: ",j,"\n")
         for (ii in 1:k){
            cat(Smtx[ii,],"\n")
         }
         if(level){
            cat("Correlations:","\n")
            print(Sj,digits=3)
         }
         ## end of if-(output) statement
      }
      wk[,(j+1)]=c(Sj)
   }
   ##
   if(output){
      par(mfcol=c(k,k))
      #### Set k0 = 4 for plotting purpose
      k0=4
      if(k > k0)par(mfcol=c(k0,k0))
      tdx=c(0,1:lags)
      jcnt=0
      if(k > 10){
         print("Skip the plots due to high dimension!")
      }
      else{
         for (j in 1:ksq){
            plot(tdx,wk[j,],type='h',xlab='lag',ylab='ccf',ylim=c(-1,1))
            abline(h=c(0))
            crit=2/sqrt(nT)
            abline(h=c(crit),lty=2)
            abline(h=c(-crit),lty=2)
            jcnt=jcnt+1
            if((jcnt==k0^2) && (k > k0)){
               jcnt=0
               cat("Hit Enter for more plots:","\n")
               readline()
            }
         }
       }
      par(mfcol=c(1,1))
      cat("Hit Enter for p-value plot of individual ccm: ","\n")
      readline()
      ## end of if-(output) statement
   }
   ## The following p-value plot was added on May 16, 2012 by Ruey Tsay.
   ### Obtain a p-value plot of ccm matrix
   r0i=solve(S0)
   R0=kronecker(r0i,r0i)
   pv=rep(0,lags)
   for (i in 1:lags){
      tmp=matrix(wk[,(i+1)],ksq,1)
      tmp1=R0%*%tmp
      ci=crossprod(tmp,tmp1)*nT*nT/(nT-i)
      pv[i]=1-pchisq(ci,ksq)
    }
   if(output){
      plot(pv,xlab='lag',ylab='p-value',ylim=c(0,1))
      abline(h=c(0))
      abline(h=c(0.05),col="blue")
      title(main="Significance plot of CCM")
    }
   ccm <- list(ccm=wk,pvalue=pv)
}


"VARorder" <- function(x,maxp=13,output=T){
   # Compute the AIC, BIC, HQ values and M-stat
   ##### Use the same number of data points in model comparison
   x1=as.matrix(x)
   nT=nrow(x1)
   k=ncol(x1)
   ksq=k*k
   if(maxp < 1)maxp=1
   enob=nT-maxp
   y=x1[(maxp+1):nT,]
   ist=maxp+1
   xmtx=cbind(rep(1,enob),x1[maxp:(nT-1),])
   if(maxp > 1){
      for (i in 2:maxp){
         xmtx=cbind(xmtx,x1[(ist-i):(nT-i),])
      }
   }
   #### in the above, y is the dependent variable, xmtx is the x-matrix
   chidet=rep(0,(maxp+1))
   s=cov(y)*(enob-1)/enob
   chidet[1]=log(det(s))
   aic=rep(0,(maxp+1))
   aic[1]=chidet[1]
   bic=aic
   hq=aic
   y=as.matrix(y)
   #
   for (p in 1:maxp){
      idm=k*p+1
      xm=xmtx[,1:idm]
      xm=as.matrix(xm)
      xpx <- crossprod(xm,xm)
      xpy <- crossprod(xm,y)
      beta <- solve(xpx,xpy)
      yhat <- xm%*%beta
      resi <- y-yhat
      sse <- crossprod(resi,resi)/enob
      #print(paste("For p = ",p,"residual variance is", sse))
      d1=log(det(sse))
      aic[p+1]=d1+(2*p*ksq)/nT
      bic[p+1]=d1+(log(nT)*p*ksq)/nT
      hq[p+1]=d1+(2*log(log(nT))*p*ksq)/nT
      chidet[p+1]=d1
   }
   maic=min(aic)
   aicor=c(1:(maxp+1))[aic==maic]-1
   mbic=min(bic)
   bicor=c(1:(maxp+1))[bic==mbic]-1
   mhq=min(hq)
   hqor=c(1:(maxp+1))[hq==mhq]-1
   
   Mstat=rep(0,maxp)
   pv=rep(0,maxp)
   for (j in 1:maxp){
      Mstat[j]=(nT-maxp-k*j-1.5)*(chidet[j]-chidet[j+1])
      pv[j]=1-pchisq(Mstat[j],ksq)
   }
   if(output){
      cat("selected order: aic = ",aicor,"\n")
      cat("selected order: bic = ",bicor,"\n")
      cat("selected order: hq = ",hqor,"\n")
      #cat("M statistic and its p-value","\n") ## comment out as
      #tmp=cbind(Mstat,pv)                     ## results shown in summary
      ##print(tmp,digits=4)
      #print(round(tmp,4))
   }
   if(output){
      n1=length(aic)-1
      ### print summary table
      cri=cbind(c(0:n1),aic,bic,hq,c(0,Mstat),c(0,pv))
      colnames(cri) <- c("p","AIC","BIC","HQ","M(p)","p-value")
      cat("Summary table: ","\n")
      ##print(cri,digits=5)
      print(round(cri,4))
   }
   
   VARorder <- list(aic=aic,aicor=aicor,bic=bic,bicor=bicor,hq=hq,hqor=hqor,Mstat=Mstat,Mpv=pv)
}
################

"VARorderI" <- function(x,maxp=13,output=T){
   # Compute the AIC, BIC, HQ values and M-stat
   ##### This is a modified version of the old program in "VARorder",
   ##### which uses the same number of data points.
   #####  This version was adopted on September 8, 2012 in Singapore.
   #####
   x1=as.matrix(x)
   nT=nrow(x1)
   k=ncol(x1)
   ksq=k*k
   if(maxp < 1)maxp=1
   ### initialization
   chidet=rep(0,(maxp+1))
   ### start with VAR(0) model, which uses just the sample means.
   s=cov(x1)*(nT-1)/nT
   chidet[1]=log(det(s))
   aic=chidet; bic=aic; hq=aic
   #
   for (p in 1:maxp){
      idm=k*p+1
      ist=p+1
      enob=nT-p
      y=as.matrix(x1[ist:nT,])
      xmtx=rep(1,enob)
      for (j in 1:p){
         xmtx=cbind(xmtx,x1[(ist-j):(nT-j),])
      }
      xm=as.matrix(xmtx)
      xpx <- crossprod(xm,xm)
      xpy <- crossprod(xm,y)
      beta <- solve(xpx,xpy)
      yhat <- xm%*%beta
      resi <- y-yhat
      sse <- crossprod(resi,resi)/enob
      #print(paste("For p = ",p,"residual variance is", sse))
      d1=log(det(sse))
      aic[p+1]=d1+(2*p*ksq)/enob
      bic[p+1]=d1+(log(enob)*p*ksq)/enob
      hq[p+1]=d1+(2*log(log(enob))*p*ksq)/enob
      chidet[p+1]=d1
   }
   maic=min(aic)
   aicor=c(1:(maxp+1))[aic==maic]-1
   mbic=min(bic)
   bicor=c(1:(maxp+1))[bic==mbic]-1
   mhq=min(hq)
   hqor=c(1:(maxp+1))[hq==mhq]-1
   
   Mstat=rep(0,maxp)
   pv=rep(0,maxp)
   for (j in 1:maxp){
      Mstat[j]=(nT-maxp-k*j-1.5)*(chidet[j]-chidet[j+1])
      pv[j]=1-pchisq(Mstat[j],ksq)
   }
   
   if(output){
      cat("selected order: aic = ",aicor,"\n")
      cat("selected order: bic = ",bicor,"\n")
      cat("selected order: hq = ",hqor,"\n")
      cat("M statistic and its p-value","\n")
      tmp=cbind(Mstat,pv)
      print(tmp,digits=4)
      ##
      n1=length(aic)-1
      ### print summary table
      cri=cbind(c(0:n1),aic,bic,hq,c(0,Mstat),c(0,pv))
      colnames(cri) <- c("p","AIC","BIC","HQ","M(p)","p-value")
      cat("Summary table: ","\n")
      print(cri,digits=5)
   }
   
   VARorderI <- list(aic=aic,aicor=aicor,bic=bic,bicor=bicor,hq=hq,hqor=hqor,Mstat=Mstat,Mpv=pv)
}
##############################


"VARpsi" <- function(Phi,lag=5){
   # Computes the psi-weight matrices of a VAR(p) model.
   # Phi=[phi1,phi2,phi3, ....] coefficient matrix
   # Created by Ruey S. Tsay, April 2009 & modified on April 2011.
   #
   # Compute MA representions
   k=nrow(Phi)
   m=ncol(Phi)
   p=floor(m/k)
   Si=diag(rep(1,k))
   if(p < 1) p =1
   if(lag < 1) lag=1
   #
   for (i in 1:lag){
      if (i <= p){
         idx=(i-1)*k
         tmp=Phi[,(idx+1):(idx+k)]
      }
      else{
         tmp=matrix(0,k,k)
      }
      #
      jj=i-1
      jp=min(jj,p)
      if(jp > 0){
         for(j in 1:jp){
            jdx=(j-1)*k
            idx=(i-j)*k
            w1=Phi[,(jdx+1):(jdx+k)]
            w2=Si[,(idx+1):(idx+k)]
            tmp=tmp+w1%*%w2
            ##print(tmp,digits=4)
         }
      }
      Si=cbind(Si,tmp)
   }
   
   VARpsi <- list(psi=Si)
}

"VARpred" <- function(model,h=1,orig=0,Out.level=F){
   # Computes the i=1, 2, ..., h-step ahead predictions of a VAR(p) model.
   # Phi=[phi1,phi2,phi3, ....] coefficient matrix
   # cnst= constant term
   # Created by Ruey S. Tsay in April 2011.
   #
   # Modified on April 20, 2011
   # It needs the program VARpsi.R to obtain the psi-weights
   # First compute the psi-weights of the VAR(p) model.
   #
   # Modifed on March 29, 2012 to include MSE of using
   #  estimated parameters.
   #
   # model is a VAR output object.
   # Out.level : control the details of output.
   #
   x=model$data
   Phi=model$Phi
   sig=model$Sigma
   Ph0=model$Ph0
   p=model$order
   cnst=model$cnst
   np=dim(Phi)[2]
   k=dim(x)[2]
   #
   nT=dim(x)[1]
   k=dim(x)[2]
   if(orig <= 0)orig=nT
   if(orig > nT)orig=nT
   psi=VARpsi(Phi,h)$psi
   beta=t(Phi)
   if(length(Ph0) < 1)Ph0=rep(0,k)
   if(p > orig){
      cat("Too few data points to produce forecasts","\n")
   }
   pred=NULL
   se=NULL
   MSE=NULL
   mse=NULL
   px=as.matrix(x[1:orig,])
   Past=px[orig,]
   if(p > 1){
      for (j in 1:(p-1)){
         Past=c(Past,px[(orig-j),])
      }
   }
   #
   # Setup to compute MSE (due to estimated parameters.
   # Compute G-matrix and construct P-matrix
   cat("orig ",orig,"\n")
   ne=orig-p
   xmtx=NULL
   P=NULL
   if(cnst)xmtx=rep(1,ne)
   xmtx=cbind(xmtx,x[p:(orig-1),])
   ist=p+1
   if(p > 1){
      for (j in 2:p){
         xmtx=cbind(xmtx,x[(ist-j):(orig-j),])
      }
   }
   xmtx=as.matrix(xmtx)
   G=t(xmtx)%*%xmtx/ne
   Ginv=solve(G)
   ##cat("G-matrix: ","\n")
   ##print(G)
   ##cat("Ginv","\n")
   ##print(Ginv)
   #
   P = Phi
   vv=Ph0
   if(p > 1){
      II=diag(rep(1,k*(p-1)))
      II=cbind(II,matrix(0,(p-1)*k,k))
      P=rbind(P,II)
      vv=c(vv,rep(0,(p-1)*k))
   }
   if(cnst){
      c1=c(1,rep(0,np))
      P=cbind(vv,P)
      P=rbind(c1,P)
   }
   ##
   ##cat("P-matrix","\n")
   ##print(P)
   #
   Sig=sig
   n1=dim(P)[2]
   MSE= (n1/orig)*sig
   for (j in 1:h){
      tmp=Ph0+matrix(Past,1,np)%*%beta
      px=rbind(px,tmp)
      if(np > k){
         Past=c(tmp,Past[1:(np-k)])
      }else{
         Past=tmp
      }
      #
      #### Compute variance of forecast errors for j > 1.
      if(j > 1){
         idx=(j-1)*k
         wk=psi[,(idx+1):(idx+k)]
         Sig=Sig+wk%*%sig%*%t(wk)
      }
      ##### Compute MSE of forecast errors for j > 1.
      if(j > 1){
         for (ii in 0:(j-1)){
            psii=diag(rep(1,k))
            if(ii > 0){
               idx=ii*k
               psii=psi[,(idx+1):(idx+k)]
            }
            P1=P^(j-1-ii)%*%Ginv
            for (jj in 0:(j-1)){
               psij=diag(rep(1,k))
               if(jj > 0){
                  jdx=jj*k
                  psij=psi[,(jdx+1):(jdx+k)]
               }
               P2=P^(j-1-jj)%*%G
               k1=sum(diag(P1%*%P2))
               MSE=(k1/orig)*psii%*%sig%*%t(psij)
            }
         }
         #
      }
      #
      se=rbind(se,sqrt(diag(Sig)))
      if(Out.level){
         cat("Covariance matrix of forecast errors at horizon: ",j,"\n")
         print(Sig)
         cat("Omega matrix at horizon: ",j,"\n")
         print(MSE)
      }
      #
      MSE=MSE+Sig
      mse=rbind(mse,sqrt(diag(MSE)))
   }
   cat("Forecasts at origin: ",orig,"\n")
   print(px[(orig+1):(orig+h),],digits=4)
   cat("Standard Errors of predictions: ","\n")
   print(se[1:h,],digits=4)
   pred=px[(orig+1):(orig+h),]
   cat("Root mean square errors of predictions: ","\n")
   print(mse[1:h,],digits=4)
   if(orig < nT){
      cat("Observations, predicted values,     errors, and MSE","\n")
      tmp=NULL
      jend=min(nT,(orig+h))
      for (t in (orig+1):jend){
         case=c(t,x[t,],px[t,],x[t,]-px[t,])
         tmp=rbind(tmp,case)
      }
      colnames(tmp) <- c("time",rep("obs",k),rep("fcst",k),rep("err",k))
      idx=c(1)
      for (j in 1:k){
         idx=c(idx,c(0,1,2)*k+j+1)
      }
      tmp = tmp[,idx]
      #print(tmp,digits=3)
      print(round(tmp,4))
   }
   
   VARpred <- list(pred=pred,se.err=se,mse=mse)
}



"VARfore" <- function(model,h=1,orig=0){
   # Computes the i=1, 2, ..., h-step ahead predictions of a VAR(p) model.
   # Phi=[phi1,phi2,phi3, ....] coefficient matrix
   # cnst= constant term
   # model is a VAR output object.
   #
   x=model$data
   Phi=model$Phi
   sig=model$Sigma
   Ph0=model$Ph0
   p=model$order
   np=dim(Phi)[2]
   #
   nT=dim(x)[1]
   k=dim(x)[2]
   if(orig <= 0)orig=nT
   if(orig > nT)orig=nT
   psi=VARpsi(Phi,h)$psi
   beta=t(Phi)
   if(length(Ph0) < 1)Ph0=rep(0,k)
   if(p > orig){
      cat("Too few data points to produce forecasts","\n")
   }
   pred=NULL
   se=NULL
   px=as.matrix(x[1:orig,])
   Past=px[orig,]
   if(p > 1){
      for (j in 1:(p-1)){
         Past=c(Past,px[(orig-j),])
      }
   }
   #
   for (j in 1:h){
      tmp=Ph0+matrix(Past,1,np)%*%beta
      px=rbind(px,tmp)
      if(np > k){
         Past=c(tmp,Past[1:(np-k)])
      }else{
         Past=tmp
      }
      Sig=sig
      if (j > 1){
         for (ii in 1:(j-1)){
            idx=ii*k
            wk=psi[,(idx+1):(idx+k)]
            Sig=Sig+wk%*%sig%*%t(wk)
         }
      }
      se=rbind(se,sqrt(diag(Sig)))
   }
   cat("Forecasts at origin: ",orig,"\n")
   print(px[(orig+1):(orig+h),],digits=4)
   cat("Standard Errors of predictions: ","\n")
   print(se[1:h,],digits=4)
   pred=px[(orig+1):(orig+h),]
   if(orig < nT){
      cat("Observations, predicted values, and errors","\n")
      tmp=NULL
      jend=min(nT,(orig+h))
      for (t in (orig+1):jend){
         case=c(t,x[t,],px[t,],x[t,]-px[t,])
         tmp=rbind(tmp,case)
      }
      colnames(tmp) <- c("time",rep("obs",k),rep("fcst",k),rep("err",k))
      idx=c(1)
      for (j in 1:k){
         idx=c(idx,c(0,1,2)*k+j+1)
      }
      tmp = tmp[,idx]
      ##print(tmp,digits=3)
      print(round(tmp,4))
   }
   VARfore <- list(pred=pred,se.err=se)
 }

"VARMAsim" <- function(nobs,arlags=NULL,malags=NULL,cnst=NULL,phi=NULL,theta=NULL,skip=200,sigma){
   # Generate VARMA(p,q) time series using Gaussian innovations.
   # p: ar order (lags can be skipped)
   # q: ma order (lags can be skipped)
   # nobs: sample size
   # cnst: constant vector
   # phi: store AR coefficient matrices [phi1,phi2,...]
   # theta: store MA coefficient matrices [theta1,theta2,...]
   # arlags: order for each AR coefficient matrix
   # malags: order for each MA coefficient matrix.
   #
   if(!is.matrix(sigma))sigma=as.matrix(sigma)
   k=nrow(sigma)
   nT=nobs+skip
   at=rmvnorm(nT,rep(0,k),sigma)
   nar=length(arlags)
   p=0
   if(nar > 0){
      arlags=sort(arlags)
      p=arlags[nar]
   }
   q=0
   nma=length(malags)
   if(nma > 0){
      malags=sort(malags)
      q=malags[nma]
   }
   ist=max(p,q)+1
   zt=matrix(0,nT,k)
   if(length(cnst)==0)cnst=rep(0,k)
   
   for (it in ist:nT){
      tmp=matrix(at[it,],1,k)
      if(nma > 0){
         for (j in 1:nma){
            jdx=(j-1)*k
            thej=theta[,(jdx+1):(jdx+k)]
            atm=matrix(at[it-malags[j],],1,k)
            tmp=tmp-atm%*%t(thej)
         }
      }
      if(nar > 0){
         for (i in 1:nar){
            idx=(i-1)*k
            phj = phi[,(idx+1):(idx+k)]
            ztm=matrix(zt[it-arlags[i],],1,k)
            tmp=tmp+ztm%*%t(phj)
         }
      }
      zt[it,]=cnst+tmp
   }
   # skip the first "skip" points
   zt=zt[(1+skip):nT,]
   at=at[(1+skip):nT,]
   
   VARMAsim <- list(series=zt,noises=at)
}

"VARirf" <- function(Phi,Sig,lag=12,orth=TRUE){
   # Computes impulse response function of a given VAR(p) model.
   # Phi: k by kp matrix of AR ceofficients, i.e. [AR1,AR2,AR3, ..., ARp]
   # Sig: residual covariance matrix
   # Output: (a) Plot and (b) Impulse response function [Psi1,Psi2, ....]
   if(!is.matrix(Phi))Phi=as.matrix(Phi)
   if(!is.matrix(Sig))Sig=as.matrix(Sig)
   # Compute MA representions: This gives impulse response function without considering Sigma.
   k=nrow(Phi)
   m=ncol(Phi)
   p=floor(m/k)
   Si=diag(rep(1,k))
   wk=c(Si)
   ## acuwk: accumulated response
   awk=c(wk)
   acuwk=c(awk)
   if(p < 1) p =1
   if(lag < 1) lag=1
   #
   for (i in 1:lag){
      if (i <= p){
         idx=(i-1)*k
         tmp=Phi[,(idx+1):(idx+k)]
      }
      else{
         tmp=matrix(0,k,k)
      }
      #
      jj=i-1
      jp=min(jj,p)
      if(jp > 0){
         for(j in 1:jp){
            jdx=(j-1)*k
            idx=(i-j)*k
            w1=Phi[,(jdx+1):(jdx+k)]
            w2=Si[,(idx+1):(idx+k)]
            tmp=tmp+w1%*%w2
            ##print(tmp,digits=4)
         }
      }
      Si=cbind(Si,tmp)
      wk=cbind(wk,c(tmp))
      awk=awk+c(tmp)
      acuwk=cbind(acuwk,awk)
      ##print(Si,digits=3)
   }
   # Compute the impulse response of orthogonal innovations
   orSi=NULL
   wk1=NULL
   awk1=NULL
   acuwk1=NULL
   if(orth){
      m1=chol(Sig)
      P=t(m1)
      wk1=cbind(wk1,c(P))
      awk1=wk1
      acuwk1=wk1
      orSi=cbind(orSi,P)
      for(i in 1:lag){
         idx=i*k
         w1=Si[,(idx+1):(idx+k)]
         w2=w1%*%P
         orSi=cbind(orSi,w2)
         wk1=cbind(wk1,c(w2))
         awk1=awk1+c(w2)
         acuwk1=cbind(acuwk1,awk1)
      }
   }
   
   tdx=c(1:(lag+1))-1
   par(mfcol=c(k,k),mai=c(0.3,0.3,0.3,0.3))
   if(orth){
      gmax=max(wk1)
      gmin=min(wk1)
      cx=(gmax-gmin)/10
      gmax=gmax+cx
      gmin=gmin-cx
      for (j in 1:k^2){
         plot(tdx,wk1[j,],type='l',xlab='lag',ylab='IRF',ylim=c(gmin,gmax),cex.axis=0.8)
         points(tdx,wk1[j,],pch='*',cex=0.8)
         title(main='Orth. innovations')
      }
      cat("Press return to continue ","\n")
      readline()
      gmax=max(acuwk1)
      gmin=min(acuwk1)
      cx=(gmax-gmin)/10
      gmax=gmax+cx
      gmin=gmin-cx
      for (j in 1:k^2){
         plot(tdx,acuwk1[j,],type='l',xlab='lag',ylab="Acu-IRF",ylim=c(gmin,gmax),cex.axis=0.8)
         points(tdx,acuwk1[j,],pch="*",cex=0.8)
         title(main='Orth. innovations')
      }
   }
   else{
      gmax=max(wk)
      gmin=min(wk)
      cx=(gmax-gmin)/10
      gmax=gmax+cx
      gmin=gmin-cx
      for(j in 1:k^2){
         plot(tdx,wk[j,],type='l',xlab='lag',ylab='IRF',ylim=c(gmin,gmax),cex.axis=0.8)
         points(tdx,wk[j,],pch='*',cex=0.8)
         title(main="Orig. innovations")
      }
      cat("Press return to continue ","\n")
      readline()
      gmax=max(acuwk)
      gmin=min(acuwk)
      cx=(gmax-gmin)/10
      gmax=gmax+cx
      gmin=gmin-cx
      for(j in 1:k^2){
         plot(tdx,acuwk[j,],type='l',xlab='lag',ylab='Acu-IRF',ylim=c(gmin,gmax),cex.axis=0.8)
         points(tdx,acuwk[j,],pch='*',cex=0.8)
         title(main="Orig. innovations")
      }
   }
   
   VARirf <- list(irf=Si,orthirf=orSi)
}

"mq" <- function(x,lag=24,adj=0){
   # Compute multivariate Ljung-Box test statistics
   #
   # adj: adjustment for the degrees of freedomm in the chi-square distribution.
   # adj is the number of coefficient parameters used in the fitted model, if any.
   #
   if(!is.matrix(x))x=as.matrix(x)
   nr=nrow(x)
   nc=ncol(x)
   g0=var(x)
   ginv=solve(g0)
   qm=0.0
   QM=NULL
   df = 0
   for (i in 1:lag){
      x1=x[(i+1):nr,]
      x2=x[1:(nr-i),]
      g = cov(x1,x2)
      g = g*(nr-i-1)/(nr-1)
      h=t(g)%*%ginv%*%g%*%ginv
      qm=qm+nr*nr*sum(diag(h))/(nr-i)
      df=df+nc*nc
      dff= df-adj
      mindeg=nc^2-1
      pv = 1
      if(dff > mindeg)pv=1-pchisq(qm,dff)
      QM=rbind(QM,c(i,qm,dff,pv))
   }
   pvs=QM[,4]
   dimnames(QM) = list(names(pvs),c("  m  ","    Q(m) ","   df  "," p-value"))
   cat("Ljung-Box Statistics: ","\n")
   printCoefmat(QM,digits = 3)
   #
   par(mfcol=c(1,1))
   plot(pvs,ylim=c(0,1),xlab="m",ylab="prob",main="p-values of Ljung-Box statistics")
   abline(h=c(0))
   lines(rep(0.05,lag),lty=2,col='blue')
}

###
"VMAorder" <- function(x,lag=20){
   # Compute multivariate Ljung-Box test statistics
   # to identify the VMA order.
   #
   if(!is.matrix(x))x=as.matrix(x)
   nr=dim(x)[1]
   nc=dim(x)[2]
   g0=var(x)
   ginv=solve(g0)
   qm=NULL
   for (i in 1:lag){
      x1=x[(i+1):nr,]
      x2=x[1:(nr-i),]
      g = cov(x1,x2)
      g = g*(nr-i-1)/(nr-1)
      h=t(g)%*%ginv%*%g%*%ginv
      qmi=nr*nr*sum(diag(h))/(nr-i)
      qm=c(qmi,qm)
   }
   tst=rev(cumsum(qm))
   ksq=nc*nc; df=ksq*lag
   QM=NULL
   for (i in 1:lag){
      pv=1-pchisq(tst[i],df)
      QM=rbind(QM,c(i,tst[i],pv))
      df=df-ksq
   }
   pvs=QM[,3]
   dimnames(QM) = list(names(pvs),c("  j  ","  Q(j,m) "," p-value"))
   cat("Q(j,m) Statistics: ","\n")
   printCoefmat(QM,digits = 3)
   ## Plot
   par(mfcol=c(1,1))
   plot(pvs,ylim=c(0,1),xlab='j',ylab='prob',main="p-values: Q(j,m) Statistics")
   abline(h=c(0))
   lines(rep(0.05,lag),lty=2,col='blue')
   #
}


#### VMA programs
"VMA" <- function(da,q=1,include.mean=T,fixed=NULL,beta=NULL,sebeta=NULL,prelim=F,details=F,thres=2.0){
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

 LLKvma <- function(par,zt=zt,q=q,fixed=fix1,include.mean=include.mean){
   ## the model used is
   ## x_t' = mu' + a_t' -a_{t-1}'theta_1' - a_{t-2}'theta_2' - ...
   ## a_t' = x_t' - mu' + a_{t-1}'theta_1'+a_{t-2}'theta_2' + ....
   k=ncol(zt)
   nT=nrow(zt)
   mu=rep(0,k)
   icnt=0; VMAcnt <- 0
   fix <- fixed
   #
   iist=0
   if(include.mean){
      iist=1
      jdx=c(1:k)[fix[1,]==1]
      icnt=length(jdx); VMAcnt <- icnt
      if(icnt > 0)
      mu[jdx]=par[1:icnt]
   }
   ### remove the mean
   for (j in 1:k){
      zt[,j]=zt[,j]-mu[j]
   }
   ## recursively compute the residual series: at
   kq=k*q
   Theta=matrix(0,kq,k)
   for (j in 1:k){
      idx=c(1:kq)[fix[(iist+1):(iist+kq),j]==1]
      jcnt=length(idx)
      if(jcnt > 0){
         Theta[idx,j]=par[(icnt+1):(icnt+jcnt)]
         icnt=icnt+jcnt
      }
   }
   # Theta = rbind[theta_1',theta_2', ..., theta_q']
   ### Checking the invertibility of t(Theta)
   TH=t(Theta)
   if(q > 1){
      tmp=cbind(diag(rep(1,(q-1)*k)),matrix(0,(q-1)*k,k))
      TH=rbind(TH,tmp)
   }
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
      ist=0
      if(VMAcnt > 0)ist=1
      for (j in 1:k){
         idx=c(1:kq)[fix[(ist+1):(ist+kq),j]==1]
         jcnt=length(idx)
         if(jcnt > 0){
            par[(icnt+1):(icnt+jcnt)]=TH[j,idx]
            icnt=icnt+jcnt
         }
      }
      ##
   }
   ##
   at=mFilter(zt,t(Theta))
   #
   sig=t(at)%*%at/nT
   ##ll=dmnorm(at,rep(0,k),sig)
   ll=dmvnorm(at,rep(0,k),sig)
   LLKvma=-sum(log(ll))
   LLKvma
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
      fit = nlminb(start = ParMA, objective = LLKvma,zt=da,fixed=fixed,include.mean=include.mean,q=q,
      lower = lowerBounds, upper = upperBounds, control = list(trace=3))
   }
   else {
      fit = nlminb(start = ParMA, objective = LLKvma, zt=da, fixed=fixed, include.mean=include.mean, q=q, 
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
         Hessian[i, j] = (LLKvma(x1,zt=da,q=q,fixed=fixed,include.mean=include.mean)
                          -LLKvma(x2,zt=da,q=q,fixed=fixed,include.mean=include.mean)
                          -LLKvma(x3,zt=da,q=q,fixed=fixed,include.mean=include.mean)
                          +LLKvma(x4,zt=da,q=q,fixed=fixed,include.mean=include.mean))/
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
   
   VMA <- list(data=da,MAorder=q,cnst=include.mean,coef=TH,secoef=seTH,residuals=at,Sigma=sig,Theta=Theta,mu=cnt,aic=aic,bic=bic)
}


####
"refVMA" <- function(model,thres=1.0){
   # This program refines the fitted models of VMA output by removing
   # insigificant parameters with abs(t-ratio) < thres.
   # model: output object from VMA
   # thres: threshold value
   #
   x = model$data
   q = model$MAorder
   cnst = model$cnst
   coef=as.matrix(model$coef)
   secoef=as.matrix(model$secoef)
   nr=dim(coef)[1]
   nc=dim(coef)[2]
   for (i in 1:nc){
      idx=is.na(secoef[,i])
      jdx=c(1:nr)[idx==T]
      secoef[jdx,i]=0.01
   }
   fix=matrix(0,nr,nc)
   for (j in 1:nc){
      tt=coef[,j]/secoef[,j]
      idx=c(1:nr)[abs(tt) >= thres]
      fix[idx,j]=1
   }
   ### Try to keep the constant if the t-ratio is greater than 1.
   if(cnst){
      tt=coef[1,]/secoef[1,]
      idx=c(1:nc)[abs(tt) > 1.0]
      if(length(idx) > 0)fix[1,idx]=1
   }
     
   mm=VMA(x,q=q,include.mean=cnst,fixed=fix,beta=coef,sebeta=secoef)
   
   refVMA <- list(data=x,MAorder=q,cnst=cnst,coef=mm$coef,secoef=mm$secoef,residuals=mm$residuals,Sigma=mm$Sigma,aic=mm$aic,bic=mm$bic,mu=mm$mu,Theta=mm$Theta)
  }

####
"VMAs" <- function(da,malags,include.mean=T,fixed=NULL,prelim=F,details=F,thres=2.0){
   # Estimation of a vector MA model using conditional MLE (Gaussian dist)
   # The MA lags are given specifically.
   #
   if(!is.matrix(da))da=as.matrix(da)
   nT <- dim(da)[1]; k <- dim(da)[2]
   nlags=length(malags)
   if(nlags < 1){
      malags=c(1)
      nlags=1
     }
   MAlag <- sort(malags)
   kq=k*nlags
   # find the maximum MA order
   q=MAlag[nlags]
 #
 THinis <- function(y,x,MAlag,include.mean){
   # use residuals of a long VAR model to obtain initial estimates of
   # VMA coefficients.
   if(!is.matrix(y))y=as.matrix(y)
   if(!is.matrix(x))x=as.matrix(x)
   nT=dim(y)[1]
   k=dim(y)[2]
   nlags=length(MAlag)
   q=MAlag[nlags]
   ist=1+q
   ne=nT-q
   if(include.mean){
      xmtx=matrix(1,ne,1)
   }
   else {
      xmtx=NULL
   }
   ymtx=y[ist:nT,]
   for (j in 1:nlags){
      jj=MAlag[j]
      xmtx=cbind(xmtx,x[(ist-jj):(nT-jj),])
   }
   xtx=crossprod(xmtx,xmtx)
   xty=crossprod(xmtx,ymtx)
   xtxinv=solve(xtx)
   beta=xtxinv%*%xty
   # compute standard errors
   resi=ymtx - xmtx%*%beta
   sse=crossprod(resi,resi)/ne
   dd=diag(xtxinv)
   sebeta=NULL
   for (j in 1:k){
      se=sqrt(dd*sse[j,j])
      sebeta=cbind(sebeta,se)
     }   
   THinis <- list(estimates=beta, se=sebeta)
  }

   # Obtain initial parameter estimates
   ### Use VAR approximation to obtain initial parameter estimates
   m1=VARorder(da,q+10,output=FALSE)
   porder=m1$aicor
   m2=VAR(da,porder,output=FALSE)
   y=da[(porder+1):nT,]
   x=m2$residuals
   m3=THinis(y,x,MAlag,include.mean)
   beta=m3$estimates
   sebeta=m3$se
   nr=dim(beta)[1]
   #
   if(prelim){
      fixed=matrix(0,nr,k)
      for (j in 1:k){
         tt=beta[,j]/sebeta[,j]
         idx=c(1:nr)[abs(tt) >= thres]
         fixed[idx,j]=1
      }
   }
   ####
   if(length(fixed)==0){fixed=matrix(1,nr,k)}
   #
   par=NULL
   separ=NULL
   ist=0
   if(include.mean){
      jdx=c(1:k)[fixed[1,]==1]
      if(length(jdx) > 0){
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
   #####
   for (j in 1:k){
      idx=c(1:(nr-ist))[fixed[(ist+1):nr,j]==1]
      if(length(idx)>0){
         par=c(par,TH[idx,j])
         separ=c(separ,seTH[idx,j])
      }
   }
   cat("Initial estimates: ",round(par,4),"\n")
   ### Set up lower and upper bounds
   lowerBounds=par; upperBounds=par
   for (j in 1:length(par)){
      lowerBounds[j] = par[j]-2*separ[j]
      upperBounds[j] = par[j]+2*separ[j]
   }
   cat("Par. lower-bounds: ",round(lowerBounds,4),"\n")
   cat("Par. upper-bounds: ",round(upperBounds,4),"\n")
###
LLKvmas <- function(par,zt=da, include.mean=include.mean, MAlag=MAlag, fixed=fixed){
   ## the model used is
   ## x_t' = mu' + a_t' -a_{t-1}'theta_1' - a_{t-2}'theta_2' - ...
   ## a_t' = x_t' - mu' + a_{t-1}'theta_1'+a_{t-2}'theta_2' + ....
   k=ncol(zt)
   nT=nrow(zt)
   nlags=length(MAlag)
   q=MAlag[nlags]
   fix <- fixed 
   #
   mu=rep(0,k)
   ist=0
   icnt=0; VMAcnt <- 0
   if(include.mean){
      ist=1
      jdx=c(1:k)[fix[1,]==1]
      icnt=length(jdx); VMAcnt <- icnt
      mu[jdx]=par[1:icnt]
      ### remove the mean
      for (j in 1:k){
         zt[,j]=zt[,j]-mu[j]
      }
   }
   ## recursively compute the residual series: at
   kq=k*nlags
   Theta=matrix(0,kq,k)
   for (j in 1:k){
      idx=c(1:kq)[fix[(ist+1):(ist+kq),j]==1]
      jcnt=length(idx)
      if(jcnt > 0){
         Theta[idx,j]=par[(icnt+1):(icnt+jcnt)]
         icnt=icnt+jcnt
      }
   }
   
   # Theta = rbind[theta_1',theta_2', ..., theta_q']
   at=zt[1:MAlag[1],]
   if(MAlag[1]==1)at=matrix(at,1,k)
   if(q >= (MAlag[1]+1)){
      for(t in (MAlag[1]+1):q){
         Past=NULL
         for (ii in 1:nlags){
            jj=MAlag[ii]
            if((t-jj) > 0){
               Past=c(Past,at[t-jj,])
            }
            else {
               Past=c(Past,rep(0,k))
            }
         }
         tmp=zt[t,]+matrix(Past,1,kq)%*%Theta
         at=rbind(at,tmp)
      }
      #end of if(q >= (MAlag[1]+1)) statement
   }
   for(t in (q+1):nT){
      Past=NULL
      for (ii in 1:nlags){
         jj=MAlag[ii]
         Past=c(Past,at[t-jj,])
      }
      tmp=zt[t,]+matrix(Past,1,kq)%*%Theta
      at=rbind(at,tmp)
   }
   #
   sig=t(at)%*%at/nT
   ##ll=dmnorm(at,rep(0,k),sig)
   ll=dmvnorm(at,rep(0,k),sig)
   LLKvmas=-sum(log(ll))
   LLKvmas
  }

   # Step 5: Estimate Parameters and Compute Numerically Hessian:
   if(details){
      fit = nlminb(start = par, objective = LLKvmas, zt=da, include.mean=include.mean, MAlag=MAlag, fixed=fixed, 
      lower = lowerBounds, upper = upperBounds, control = list(trace=3))
   }
   else {
      fit = nlminb(start = par, objective = LLKvmas, zt=da, include.mean=include.mean, MAlag=MAlag, fixed=fixed, 
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
         Hessian[i, j] = (LLKvmas(x1,zt=da,include.mean=include.mean,MAlag=MAlag,fixed=fixed)
                         -LLKvmas(x2,zt=da,include.mean=include.mean,MAlag=MAlag,fixed=fixed)
                         -LLKvmas(x3,zt=da,include.mean=include.mean,MAlag=MAlag,fixed=fixed)
                         +LLKvmas(x4,zt=da,include.mean=include.mean,MAlag=MAlag,fixed=fixed))/
         (4*epsilon[i]*epsilon[j])
      }
   }
   est=fit$par
   cat("Final    Estimates: ",est,"\n")
   # Step 6: Create and Print Summary Report:
   se.coef = sqrt(diag(solve(Hessian)))
   tval = fit$par/se.coef
   matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
   dimnames(matcoef) = list(names(tval), c(" Estimate",
   " Std. Error", " t value", "Pr(>|t|)"))
   cat("\nCoefficient(s):\n")
   printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
   
   cat("---","\n")
   cat("Estimates in matrix form:","\n")
   icnt=0
   ist=0
   cnt=rep(0,k)
   secnt=rep(1,k)
   # handle the mean,if any.
   if(include.mean){
      ist=1
      jdx=c(1:k)[fixed[1,]==1]
      icnt=length(jdx)
      if(icnt > 0){
         cnt[jdx]=est[1:icnt]
         secnt=se.coef[1:icnt]
         cat("Constant term: ","\n")
         cat("Estimates: ",cnt,"\n")
      }
   }
   cat("MA coefficient matrix","\n")
   TH=matrix(0,kq,k)
   seTH=matrix(1,kq,k)
   for (j in 1:k){
      idx=c(1:kq)[fixed[(ist+1):nr,j]==1]
      jcnt=length(idx)
      if(jcnt > 0){
         TH[idx,j]=est[(icnt+1):(icnt+jcnt)]
         seTH[idx,j]=se.coef[(icnt+1):(icnt+jcnt)]
         icnt=icnt+jcnt
      }
   }
   #
   ####print(TH)
   #
   icnt=0
   for (i in 1:nlags){
      ii=MAlag[i]
      cat("MA(",ii,")-matrix","\n")
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
   at=zt[1:MAlag[1],]
   if(MAlag[1]==1)at=matrix(at,1,k)
   if(q >= (MAlag[1]+1)){
      for(t in (MAlag[1]+1):q){
         Past=NULL
         for (ii in 1:nlags){
            jj=MAlag[ii]
            if((t-jj) > 0){
               Past=c(Past,at[t-jj,])
            }
            else {
               Past=c(Past,rep(0,k))
            }
         }
         tmp=zt[t,]+matrix(Past,1,kq)%*%TH
         at=rbind(at,tmp)
      }
      #end of the statement if(q > (MAlag[1]+1)
   }
   for(t in (q+1):nT){
      Past=NULL
      for (ii in 1:nlags){
         jj=MAlag[ii]
         Past=c(Past,at[t-jj,])
      }
      tmp=zt[t,]+matrix(Past,1,kq)%*%TH
      at=rbind(at,tmp)
   }
   #
   sig=t(at)%*%at/nT
   cat(" ","\n")
   cat("Residuals cov-matrix:","\n")
   print(sig)
   dd=det(sig)
   d1=log(dd)
   aic=d1+2*npar/nT
   bic=d1+log(nT)*npar/nT
   cat("---","\n")
   cat("aic = ",aic,"\n")
   cat("bic = ",bic,"\n")
   ###
   Theta=t(TH)
   if(include.mean){
      TH=rbind(cnt,TH)
      seTH=rbind(secnt,seTH)
   }
   
   VMAs <- list(data=da,MAlags=MAlag,cnst=include.mean,coef=TH,secoef=seTH,residuals=at,aic=aic,bic=bic,Sigma=sig,Theta=Theta,mu=cnt,MAorder=q)
}


####
"refVMAs" <- function(model,thres=2){
   # This program refines a fittd VMAs model by removing insignificant parameters defined as
   # abs(t-ratio) < thres.
   #
   # model: an output object from VMAs program
   # thres: threshold
   #
   x = model$data
   malags = model$MAlags
   cnst = model$cnst
   coef=model$coef
   secoef=model$secoef
   nr=dim(coef)[1]
   nc=dim(coef)[2]
   for (i in 1:nr){
      for (j in 1:nc){
         if(secoef[i,j] < 10^(-8))secoef[i,j]=1.0
      }
   }
   k=dim(x)[2]
   nr=dim(coef)[1]
   nc=dim(coef)[2]
   fix=matrix(0,nr,k)
   for (j in 1:k){
      tt=coef[,j]/secoef[,j]
      idx=c(1:nr)[abs(tt) >= thres]
      if(length(idx) > 0)fix[idx,j]=1
   }
   ### Try to keep the constant if the t-ratio is greater then 1.
   if(cnst){
      tt=coef[1,]/secoef[1,]
      idx=c(1:nc)[abs(tt) > 1.0]
      if(length(idx) > 0)fix[1,idx]=1
   }
   
   mm=VMAs(x,malags,include.mean=cnst,fixed=fix)
   
   refVMAs <- list(data=x,MAlags=malags,cnst=cnst,coef=mm$coef,secoef=mm$secoef,residuals=mm$residuals,aic=mm$aic,bic=mm$bic,Sigma=mm$Sigma,Theta=mm$Theta,mu=mm$mu,MAorder=mm$MAorder)
}

#####
"VMApred" <- function(model,h=1,orig=0){
   # Computes the i=1, 2, ..., h-step ahead predictions of a VMA(q) model.
   #
   # model is a VMA output object.
   # created on April 20, 2011
   #
   x=model$data
   resi=model$residuals
   Theta=model$Theta
   sig=model$Sigma
   mu=model$mu
   q=model$MAorder
   np=dim(Theta)[2]
   psi=-Theta
   #
   nT=dim(x)[1]
   k=dim(x)[2]
   if(orig <= 0)orig=nT
   if(orig > T)orig=nT
   if(length(mu) < 1)mu=rep(0,k)
   if(q > orig){
      cat("Too few data points to produce forecasts","\n")
   }
   pred=NULL
   se=NULL
   px=as.matrix(x[1:orig,])
   for (j in 1:h){
      fcst=mu
      t=orig+j
      for (i in 1:q){
         jdx=(i-1)*k
         t1=t-i
         if(t1 <= orig){
            theta=Theta[,(jdx+1):(jdx+k)]
            fcst=fcst-matrix(resi[t1,],1,k)%*%t(theta)
         }
      }
      px=rbind(px,fcst)
      #
      Sig=sig
      if (j > 1){
         jj=min(q,(j-1))
         for (ii in 1:jj){
            idx=(ii-1)*k
            wk=psi[,(idx+1):(idx+k)]
            Sig=Sig+wk%*%sig%*%t(wk)
         }
      }
      se=rbind(se,sqrt(diag(Sig)))
   }
   cat("Forecasts at origin: ",orig,"\n")
   print(px[(orig+1):(orig+h),],digits=4)
   cat("Standard Errors of predictions: ","\n")
   print(se[1:h,],digits=4)
   pred=px[(orig+1):(orig+h),]
   if(orig < nT){
      cat("Observations, predicted values, and errors","\n")
      tmp=NULL
      jend=min(nT,(orig+h))
      for (t in (orig+1):jend){
         case=c(t,x[t,],px[t,],x[t,]-px[t,])
         tmp=rbind(tmp,case)
      }
      colnames(tmp) <- c("time",rep("obs",k),rep("fcst",k),rep("err",k))
      idx=c(1)
      for (j in 1:k){
         idx=c(idx,c(0,1,2)*k+j+1)
      }
      tmp = tmp[,idx]
      print(tmp,digits=3)
   }
   
   VMApred <- list(pred=pred,se.err=se)
}

####
"VARMA" <- function(da,p=0,q=0,include.mean=T,fixed=NULL,beta=NULL,sebeta=NULL,prelim=F,details=F,thres=2.0){
   # Estimation of a vector ARMA model using conditional MLE (Gaussian dist)
   #
   # When prelim=TRUE, fixed is assigned based on the results of AR approximation.
   # Here "thres" is used only when prelim = TRUE.
   #
   if(!is.matrix(da))da=as.matrix(da)
   nT=dim(da)[1];   k=dim(da)[2]
   # basic setup.
   if(p < 0)p=0
   if(q < 0)q=0
   if((p+q) < 1)p=1
   pqmax=max(p,q)
   kq=k*q
   kp=k*p
#
 iniEST <- function(y,x,p,q,include.mean){
   # use residuals of a long VAR model to obtain initial estimates of
   # VARMA coefficients.
   if(!is.matrix(y))y=as.matrix(y)
   if(!is.matrix(x))x=as.matrix(x)
   nT=dim(y)[1]
   k=dim(y)[2]
   pq=max(p,q)
   ist=1+pq
   ne=nT-pq
   if(include.mean){
      xmtx=matrix(1,ne,1)
     }
     else {
      xmtx=NULL
     }
   ymtx=as.matrix(y[ist:nT,])
   if(p > 0){
      for (j in 1:p){
         xmtx=cbind(xmtx,y[(ist-j):(nT-j),])
       }
    }
   if(q > 0){
      for (j in 1:q){
         xmtx=cbind(xmtx,x[(ist-j):(nT-j),])
      }
     }
   xmtx=as.matrix(xmtx)
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
   iniEST <- list(estimates=beta,se=sebeta)
  }

  if(length(fixed) < 1){ 
   m1=VARorder(da,p+q+9,output=FALSE)
   porder=m1$aicor
   if(porder < 1)porder=1
   m2=VAR(da,porder,output=FALSE)
   y=da[(porder+1):nT,]
   x=m2$residuals
   m3=iniEST(y,x,p,q,include.mean)
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
    if(length(fixed)==0){fixed=matrix(1,nr,k)}
   }
   # Identify parameters to be estimated.
   par=NULL
   separ=NULL
   ist=0
   if(include.mean){
      jdx=c(1:k)[fixed[1,]==1]
      if(length(jdx) > 0){
         par=beta[1,jdx]
         separ=sebeta[1,jdx]
      }
      ist=1
   }
   if(p > 0){
      for (j in 1:k){
         idx=c(1:kp)[fixed[(ist+1):(ist+kp),j]==1]
         if(length(idx) > 0){
            tmp=beta[(ist+1):(ist+kp),j]
            setmp=sebeta[(ist+1):(ist+kp),j]
            par=c(par,tmp[idx])
            separ=c(separ,setmp[idx])
         }
         #end of j-loop
      }
      ist=ist+kp
   }
   #
   if(q > 0){
      for (j in 1:k){
         idx=c(1:kq)[fixed[(ist+1):(ist+kq),j]==1]
         if(length(idx) > 0){
            tmp=beta[(ist+1):(ist+kq),j]
            setmp=sebeta[(ist+1):(ist+kq),j]
            par=c(par,tmp[idx])
            separ=c(separ,setmp[idx])
         }
       }
   }
   #########
   cat("Number of parameters: ",length(par),"\n")
   cat("initial estimates: ",round(par,4),"\n")
   ### Set up lower and upper bounds
   lowerBounds=par; upperBounds=par
   for (j in 1:length(par)){
      lowerBounds[j] = par[j]-2*separ[j]
      upperBounds[j] = par[j]+2*separ[j]
   }
   cat("Par. lower-bounds: ",round(lowerBounds,4),"\n")
   cat("Par. upper-bounds: ",round(upperBounds,4),"\n")
#### likelihood function
LLKvarma <- function(par,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed){
   ##
   nT=dim(zt)[1]; k=dim(zt)[2]
   pqmax=max(p,q)
   kp=k*p
   kq=k*q
   ###  Assign parameters to their proper locations in the program.
   beta=NULL
   ist=0
   icnt=0
   Ph0=rep(0,k)
   if(include.mean){
      idx=c(1:k)[fixed[1,]==1]
      icnt=length(idx)
      if(icnt > 0){
         Ph0[idx]=par[1:icnt]
      }
      ist=1
      beta=rbind(beta,Ph0)
   }
   PH=NULL
   if(p > 0){
      PH = matrix(0,kp,k)
      for (j in 1:k){
         idx=c(1:kp)[fixed[(ist+1):(ist+kp),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            PH[idx,j]=par[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
         # end of j-loop
      }
      ist=ist+kp
      beta=rbind(beta,PH)
      #end of if (p > 0)
   }
   #
   TH=NULL
   if(q > 0){
      TH=matrix(0,kq,k)
      for (j in 1:k){
         idx=c(1:kq)[fixed[(ist+1):(ist+kq),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            TH[idx,j]=par[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
         # end of j-loop
      }
      beta=rbind(beta,TH)
      # end of if(q > 0).
   }
   ### recursively compute the residuals
   istart=pqmax+1
   #### consider the case t from 1 to pqmatx
   at=matrix((zt[1,]-Ph0),1,k)
   if(pqmax > 1){
      for (t in 2:pqmax){
         tmp=matrix((zt[t,]-Ph0),1,k)
         if(p > 0){
            for (j in 1:p){
               if((t-j) > 0){
                  jdx=(j-1)*k
                  tmp1=matrix(zt[(t-j),],1,k)%*%as.matrix(PH[(jdx+1):(jdx+k),])
                  tmp=tmp-tmp1
               }
               # end of j-loop
            }
            # end of if(p > 0) statement
         }
         #
         if(q > 0){
            for (j in 1:q){
               jdx=(j-1)*k
               if((t-j)>0){
                  tmp2=matrix(at[(t-j),],1,k)%*%as.matrix(TH[(jdx+1):(jdx+k),])
                  tmp=tmp-tmp2
               }
               #end of j-loop
            }
            #end of if(q > 0) statement
         }
         at=rbind(at,tmp)
      }
      # end of if(pqmax > 1) statement
   }
   ### for t from ist on
   Pcnt = NULL
   idim=kp+kq
   if(include.mean){
      Pcnt=c(1)
      idim=idim+1
   }
   for (t in istart:nT){
      Past=NULL
      if(p > 0){
         for (j in 1:p){
            Past=c(Past,zt[(t-j),])
         }
      }
      if(q > 0){
         for (j in 1:q){
            Past=c(Past,at[(t-j),])
         }
      }
      tmp = matrix(c(Pcnt,Past),1,idim)%*%beta
      tmp3=zt[t,]-tmp
      at=rbind(at,tmp3)
   }
   #### skip the first max(p,q) residuals.
   at=at[(istart:nT),]
   sig=t(at)%*%at/(nT-pqmax)
   #ll=dmnorm(at,rep(0,k),sig)
   ll=dmvnorm(at,rep(0,k),sig)
   LLKvarma=-sum(log(ll))
   LLKvarma
  }

   # Step 5: Estimate Parameters and Compute Numerically Hessian:
   if(details){
      fit = nlminb(start = par, objective = LLKvarma,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed,
      lower = lowerBounds, upper = upperBounds, control = list(trace=3))
   }
   else {
      fit = nlminb(start = par, objective = LLKvarma, zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed, 
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
         Hessian[i, j] = (LLKvarma(x1,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed)
                         -LLKvarma(x2,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed)
                         -LLKvarma(x3,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed)
                         +LLKvarma(x4,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed))/
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
   ### restore estimates to the format of unconstrained case for printing purpose.
   #### icnt: parameter count
   #### ist: location count
   ist=0
   icnt = 0
   Ph0=rep(0,k)
   sePh0=rep(0,k)
   beta=NULL
   sebeta=NULL
   if(include.mean){
      idx=c(1:k)[fixed[1,]==1]
      icnt=length(idx)
      if(icnt > 0){
         Ph0[idx]=est[1:icnt]
         sePh0[idx]=se.coef[1:icnt]
      }
      ist=1
      beta=rbind(beta,Ph0)
      sebeta=rbind(sebeta,sePh0)
   }
   PH=NULL
   sePH=NULL
   if(p > 0){
      PH=matrix(0,kp,k)
      sePH=matrix(0,kp,k)
      for (j in 1:k){
         idx=c(1:kp)[fixed[(ist+1):(ist+kp),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            PH[idx,j]=est[(icnt+1):(icnt+jdx)]
            sePH[idx,j]=se.coef[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
         # end of j-loop
      }
      #end of if (p > 0)
      ist=ist+kp
      beta=rbind(beta,PH)
      sebeta=rbind(sebeta,sePH)
   }
   #
   TH=NULL
   seTH=NULL
   if(q > 0){
      TH=matrix(0,kq,k)
      seTH=matrix(0,kq,k)
      for (j in 1:k){
         idx=c(1:kq)[fixed[(ist+1):(ist+kq),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            TH[idx,j]=est[(icnt+1):(icnt+jdx)]
            seTH[idx,j]=se.coef[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
         # end of j-loop
      }
      # end of if(q > 0).
      beta=rbind(beta,TH)
      sebeta=rbind(sebeta,seTH)
   }
   #########
   cat("---","\n")
   cat("Estimates in matrix form:","\n")
   if(include.mean){
      cat("Constant term: ","\n")
      cat("Estimates: ",Ph0,"\n")
   }
   if(p > 0){
      cat("AR coefficient matrix","\n")
      jcnt=0
      for (i in 1:p){
         cat("AR(",i,")-matrix","\n")
         ph=t(PH[(jcnt+1):(jcnt+k),])
         print(ph,digits=3)
         jcnt=jcnt+k
      }
      # end of if (p > 0)
   }
   if(q > 0){
      cat("MA coefficient matrix","\n")
      icnt=0
      for (i in 1:q){
         cat("MA(",i,")-matrix","\n")
         theta=-t(TH[(icnt+1):(icnt+k),])
         print(theta,digits=3)
         icnt=icnt+k
      }
      # end of the statement if(q > 0)
   }
   ##### Compute the residuals
   zt=da
   ist=pqmax+1
   #### consider the case t from 1 to pqmatx
   at=matrix((zt[1,]-Ph0),1,k)
   if(pqmax > 1){
      for (t in 2:pqmax){
         tmp=matrix((zt[t,]-Ph0),1,k)
         if(p > 0){
            for (j in 1:p){
               if((t-j) > 0){
                  jdx=(j-1)*k
                  tmp1=matrix(zt[(t-j),],1,k)%*%as.matrix(PH[(jdx+1):(jdx+k),])
                  tmp=tmp-tmp1
               }
               # end of j-loop
            }
            # end of if(p > 0) statement
         }
         #
         if(q > 0){
            for (j in 1:q){
               jdx=(j-1)*k
               if((t-j)>0){
                  tmp2=matrix(at[(t-j),],1,k)%*%as.matrix(TH[(jdx+1):(jdx+k),])
                  tmp=tmp-tmp2
               }
               #end of j-loop
            }
            #end of if(q > 0) statement
         }
         at=rbind(at,tmp)
         # end of for(t in 2:pqmax)
      }
      # end of if(pqmax > 1) statement
   }
   
   ### for t from ist on
   Pcnt = NULL
   idim=kp+kq
   if(include.mean){
      Pcnt=c(1)
      idim=idim+1
   }
   #
   for (t in ist:nT){
      Past=NULL
      if(p > 0){
         for (j in 1:p){
            Past=c(Past,zt[(t-j),])
         }
      }
      if(q > 0){
         for (j in 1:q){
            Past=c(Past,at[(t-j),])
         }
      }
      tmp = matrix(c(Pcnt,Past),1,idim)%*%beta
      tmp3=zt[t,]-tmp
      at=rbind(at,tmp3)
   }
   #### skip the first max(p,q) residuals.
   at=at[(ist:nT),]
   sig=t(at)%*%at/(nT-pqmax)
   ##
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
   if(length(PH) > 0)PH=t(PH)
   if(length(TH) > 0)TH=-t(TH)
   
   VARMA <- list(data=da,ARorder=p,MAorder=q,cnst=include.mean,coef=beta,secoef=sebeta,residuals=at,Sigma=sig,aic=aic,bic=bic,Phi=PH,Theta=TH,Ph0=Ph0)
}


####
"refVARMA" <- function(model,thres=1.5){
   # This program refines the fitted models of VARMA output by removing
   # insigificant parameters with abs(t-ratio) < thres.
   # model: output object from VARMA
   # thres: threshold value
   #
   x = model$data
   p1 = model$ARorder
   q1 = model$MAorder
   cnst = model$cnst
   coef=as.matrix(model$coef)
   secoef=as.matrix(model$secoef)
   nr=dim(coef)[1]
   nc=dim(coef)[2]
   for (j in 1:nc){
      idx=is.na(secoef[,j])
      jdx=c(1:nr)[idx==T]
      secoef[jdx,j]=0.01
   }
   fix=matrix(0,nr,nc)
   for (j in 1:nc){
      tt=coef[,j]/secoef[,j]
      idx=c(1:nr)[abs(tt) >= thres]
      fix[idx,j]=1
   }
   ### Try to keep the constant if the t-ratio is greater then 1.
   if(cnst){
      tt=coef[1,]/secoef[1,]
      idx=c(1:nc)[abs(tt) > 1.0]
      if(length(idx) > 0)fix[1,idx]=1
   }
   
   mm=VARMA(x,p=p1,q=q1,include.mean=cnst,fixed=fix,beta=coef,sebeta=secoef)
   
   refVARMA <- list(data=x,coef=mm$coef,secoef=mm$secoef,ARorder=p1,MAorder=q1,cnst=cnst,residuals=mm$residuals,Ph0=mm$Ph0,Phi=mm$Phi,Theta=mm$Theta,Sigma=mm$Sigma,aic=mm$aic,bic=mm$bic)
}
######
"VARMApred" <- function(model,h=1,orig=0){
   ## Compute forecasts and forecast error covariance of a VARMA mdoel.
   ## created April 21, 2011 by Ruey S. Tsay
   #
   # model: an output from VARMA command.
   #
   x=as.matrix(model$data)
   resi=as.matrix(model$residuals)
   sig=model$Sigma
   Phi=model$Phi
   Theta=model$Theta
   Ph0=model$Ph0
   p=model$ARorder
   q=model$MAorder
   #
   if(p < 0)p=0
   if(q < 0)q=0
   if(h < 1)h=1
   nT=dim(x)[1]
   k=dim(x)[2]
   T1=dim(resi)[1]
   ## In case the residuals is shorter due to conditional MLE estimation.
   if(nT > T1){
      r1=matrix(0,(nT-T1),k)
      resi=rbind(r1,resi)
   }
   #
   if(length(Ph0) < 1)Ph0=rep(0,k)
   if(orig < 1)orig=nT
   if(orig > T)orig=nT
   px=x[1:orig,]
   presi=resi[1:orig,]
   # Compute the psi-weights for the variance of forecast errors.
   psi=diag(rep(1,k))
   wk=c(psi)
   lag=max(1,h)
   #
   for (i in 1:lag){
      if (i <= p){
         idx=(i-1)*k
         tmp=Phi[,(idx+1):(idx+k)]
      }
      else{
         tmp=matrix(0,k,k)
      }
      if(i <= q){
         mdx=(i-1)*k
         tmp=tmp-Theta[,(mdx+1):(mdx+k)]
      }
      #
      jj=i-1
      jp=min(jj,p)
      if(jp > 0){
         for(j in 1:jp){
            jdx=(j-1)*k
            idx=(i-j)*k
            w1=Phi[,(jdx+1):(jdx+k)]
            w2=psi[,(idx+1):(idx+k)]
            tmp=tmp+w1%*%w2
            ##print(tmp,digits=4)
         }
      }
      psi=cbind(psi,tmp)
      wk=cbind(wk,c(tmp))
      ##print(psi,digits=3)
   }
   ### Compute the forecasts and their standard errors
   sefcst=NULL
   for (j in 1:h){
      fcst=Ph0
      Sig=sig
      t=orig+j
      ### AR part
      if(p > 0){
         for (ii in 1:p){
            idx=(ii-1)*k
            ph=Phi[,(idx+1):(idx+k)]
            fcst=fcst + matrix(px[(t-ii),],1,k)%*%t(ph)
         }
         #end of AR part
      }
      ### MA part
      if(q > 0){
         for (jj in 1:q){
            idx=(jj-1)*k
            if((t-jj) <= orig){
               th=Theta[,(idx+1):(idx+k)]
               fcst=fcst - matrix(resi[(t-jj),],1,k)%*%t(th)
            }
            # end of jj-loop
         }
         # end of MA part
      }
      px=rbind(px,fcst)
      # compute standard errors of forecasts
      if(j > 1){
         Sig=sig
         for (jj in 2:j){
            jdx=(jj-1)*k
            wk=psi[,(jdx+1):(jdx+k)]
            Sig=Sig + wk%*%sig%*%t(wk)
         }
      }
      sefcst=rbind(sefcst,sqrt(diag(Sig)))
   }
   cat("Predictions at origin ",orig,"\n")
   print(px[(orig+1):(orig+h),],digits=4)
   cat("Standard errors of predictions","\n")
   if(h == 1){
      print(sefcst,digits=4)
   }
   else {
      print(sefcst[1:h,],digits=4)
   }
   #### if orig < nT, print out actual values.
   if(orig < nT){
      cat("Observations, predictions, and errors: ","\n")
      tmp=NULL
      jend=min(nT,orig+h)
      for (t in (orig+1):jend){
         case=c(t,x[t,],px[t,],x[t,]-px[t,])
         tmp=rbind(tmp,case)
      }
      colnames(tmp) <- c("time",rep("obs",k),rep("fcst",k),rep("err",k))
      idx=c(1)
      for (j in 1:k){
         idx=c(idx,c(0,1,2)*k+j+1)
      }
      tmp = tmp[,idx]
      print(tmp,digits=4)
   }
   
   VARMApred <- list(pred=px[(orig+1):(orig+h),],se.err=sefcst,orig=orig)
   #end of the program
}


###
"VARecm" <- function(x,p=1,wt,include.const=FALSE){
   # Fits an error-correction VAR model.
   if(!is.matrix(x))x=as.matrix(x)
   nT=dim(x)[1]
   k=dim(x)[2]
   dx=x[2:nT,]-x[1:(nT-1),]
   dx=rbind(rep(0,k),dx)
   wtadj=wt-mean(wt)
   idm=k*(p-1)+1
   if(include.const)idm=idm+1
   # effective sample size
   ist=max(1,p)
   ne=nT-ist+1
   y=dx[ist:nT,]
   xmtx=wtadj[(ist-1):(nT-1)]
   if(include.const)xmtx=cbind(xmtx,rep(1,(nT-ist+1)))
   if(p > 1){
      for (i in 2:p){
         ii=i-1
         xmtx=cbind(xmtx,dx[(ist-ii):(nT-ii),])
      }
   }
   y=as.matrix(y)
   xmtx=as.matrix(xmtx)
   xpx = t(xmtx)%*%xmtx
   xpxinv=solve(xpx)
   xpy=t(xmtx)%*%y
   beta=xpxinv%*%xpy
   yhat=xmtx%*%beta
   resi=y-yhat
   sse=(t(resi)%*%resi)/ne
   alpha=beta[1,]
   icnt=1
   if(include.const){
      c=beta[2,]
      icnt=2
   }
   dd=diag(xpxinv)
   sdbeta=matrix(0,idm,k)
   for (i in 1:k){
      sdbeta[,i]=sqrt(sse[i,i]*dd)
   }
   se=sdbeta[1,]
   cat("alpha: ","\n")
   print(alpha,digits=3)
   cat("standard error","\n")
   print(se,digits=3)
   if(include.const){
      cat("constant term:","\n")
      print(c,digits=3)
      se=sdbeta[2,]
      cat("standard error","\n")
      print(se,digits=3)
   }
   
   cat("AR coefficient matrix","\n")
   jst=icnt
   for (i in 1:(p-1)){
      cat("AR(",i,")-matrix","\n")
      phi=t(beta[(jst+1):(jst+k),])
      se=t(sdbeta[(jst+1):(jst+k),])
      print(phi,digits=3)
      cat("standard error","\n")
      print(se,digits=3)
      jst=jst+k
      ###cat("      ","\n")
   }
   cat("-----","\n")
   cat("Residuals cov-mtx:","\n")
   print(sse)
   #sse=sse*ne/T
   cat("      ","\n")
   dd=det(sse)
   cat("det(sse) = ",dd,"\n")
   d1=log(dd)
   aic=d1+(2*idm*k)/nT
   bic=d1+log(nT)*idm*k/nT
   cat("AIC = ",aic,"\n")
   cat("BIC = ",bic,"\n")
   
   VARecm<-list(coef=beta,aic=aic,bic=bic,residuals=resi,secoef=sdbeta,Sigma=sse)
}

##### Mmodel checking
"MTSdiag" <- function(model,gof=24,adj=0,level=F){
   # perform model checking for a multivariate time series model.
   # m1 is a VARMA, VMA, VAR type of models.
   #
   # adj: number of coefficient parameters in the fitted model.(without counting
   #       those in the mean and covariance matrix)
   # level: switch to print residual CCM matrices.
   ###
   resi=model$residuals
   colnames(resi) <- colnames(model$data)
   ccm(resi,lags=gof,level=level)
   cat("Hit Enter to compute MQ-statistics:","\n")
   readline()
   mq(resi,lag=gof,adj=adj)
   cat("Hit Enter to obtain residual plots:","\n")
   readline()
   MTSplot(resi)
 }

####
"tfm" <- function(y,x,b=0,s=1,p=0,q=0){
   # Estimate a special transfer function model. Specifically,
   # fit an ARMA(p,q) model to [y -(w0+w1*B+w2*B**2+...+ws*B**s)B**b x].
   # b: delay
   # s: order of the transfer function polynomial.
   # note: Length(y) == length(x) & Missing values are not allowed.
   #
   # Created by R.S. Tsay, March 2009
   #
   nT=length(y)
   T1=length(x)
   nT=min(nT,T1)
   mx=b+s
   ist=mx+1
   y1=y[ist:nT]
   X=x[(s+1):(nT-b)]
   if(s > 0){
      for (i in 1:s){
         X=cbind(X,x[(s+1-i):(nT-b-i)])
      }
   }
   nx=ncol(X)
   m1=arima(y1,order=c(p,0,q),xreg=X)
   se=sqrt(diag(m1$var.coef))
   coef.arma=NULL
   se.arma=NULL
   pq=p+q
   if(pq > 0){
      coef.arma=m1$coef[1:pq]
      se.arma=se[1:pq]
      p1=cbind(coef.arma,se.arma)
      cat("ARMA coefficients & s.e.:","\n")
      print(t(p1),digits=3)
   }
   v=m1$coef[(pq+1):(pq+1+nx)]
   se.v=se[(pq+1):(pq+1+nx)]
   pr=cbind(v,se.v)
   cat("Transfer function coefficients & s.e.:","\n")
   print(t(pr),digits=3)
   res=m1$residuals
   beta=matrix(v[2:(nx+1)],nx,1)
   nt=y1-v[1]-X%*%beta
   
   tfm <- list(coef=v,se.coef=se.v,coef.arma=coef.arma,se.arma=se.arma,nt=nt,residuals=res)
  }

####
"tfm1" <- function(y,x,orderN,orderX){
   ## Estimation of a transfer function model with ONE exogenous variable
   ### The model is Y_t -c0 -w(B)/d(B)X_t = theta(B)/phi(B)a_t.
   ### orderN = c(p,d,q) for the ARMA part
   ### orderX = c(r,s,b) where d(B) = 1 - d_1B - ... - d_r B^r
   ###            and w(B) = w_0+w_1B + ... + w_s B^s and b is the delay.
   ###
   ### par=c(c0,w0,w1,...,ws,d1,...,dr,phi,theta); Feb. 2012.
   ###
   dify = orderN[2]; dY=y; dX=x
   if(dify > 0){
      dY <- y[(dify+1):length(y)]-y[1:(length(y)-dify)]
      dX <- x[(dify+1):length(x)]-x[1:(length(x)-dify)]
    }
   N = length(dY); N1=length(dX)
   if(N < N1) N1=N; if(N1 < N)N=N1
   phi=NULL; theta=NULL; ome=NULL; del=NULL
   r=orderX[1]; s=orderX[2]; b=orderX[3]; p=orderN[1]; q=orderN[3]
   r=max(r,0); s=max(0,s); b=max(0,b); p=max(0,p); q=max(0,q)
### subroutines used
 Nlike <- function(par,dY=dY,dX=dX,orderN=orderN,orderX=orderX){
   resi = Gaulike(par,dY=dY,dX=dX,orderN=orderN,orderX=orderX)
   sig=sqrt(var(resi))
   n1=length(resi)
   Nlike=-sum(log(dnorm(resi,mean=rep(0,n1),sd=sig)))
  }

 Gaulike <- function(par,dY=dY,dX=dX,orderN=orderN,orderX=orderX){
   p=orderN[1]; q=orderN[3]; r=orderX[1]; s=orderX[2]; b=orderX[3]
   c0=par[1]
   ome=par[2:(2+s)]
   #
   if(r > 0)del=par[(2+s+1):(2+s+r)]
   if(p > 0)phi=par[(3+r+s):(r+s+2+p)]
   if(q > 0)theta=par[(3+r+s+p):(2+r+s+p+q)]
   N=length(dY)
   ist=r+1
   N1t=dY-c0
   Nt = dX
   if(r > 0){
     Nt=filter(dX,del,method="r",init=rep(mean(dX),r))
     }
   ##
   ist=b+s+1
   N=length(Nt)
   N1t=N1t[ist:N]-ome[1]*Nt[(ist-b):(N-b)]
   if(s > 0){
      for (j in 1:s){
         N1t=N1t-ome[j+1]*Nt[(ist-j-b):(N-j-b)]
       }
    }
   N1=length(N1t)
   resi=N1t[(p+1):N1]
   if(p > 0){
      for (j in 1:p){
         resi=resi-phi[j]*N1t[(p+1-j):(N1-j)]
        }
     }
   #
   if(q > 0)resi=filter(resi,theta,method="r",init=rep(0,q))
   Gaulike = resi
 }

### Obtain the N(t) series
 Nts <- function(par,dY=dY,dX=dX,orderN=orderN,orderX=orderX){
   p=orderN[1]; q=orderN[3]; r=orderX[1]; s=orderX[2]; b=orderX[3]
   c0=par[1]
   ome=par[2:(2+s)]
   #
   if(r > 0)del=par[(2+s+1):(2+s+r)]
   N=length(dY)
   ist=r+1
   N1t=dY-c0
   Nt=dX
   if(r > 0){
     Nt=filter(dX,del,method="r",init=rep(mean(dX),r))
    }
   ##
   ist=b+s+1
   N=length(Nt)
   N1t=N1t[ist:N]-ome[1]*Nt[(ist-b):(N-b)]
   if(s > 0){
      for (j in 1:s){
         N1t=N1t-ome[j+1]*Nt[(ist-j-b):(N-j-b)]
        }
     }
   Nts=N1t
  }
####
   ## r = 0, the model can be fitted by the regular "arima" command.
   if(r==0){
      nobe=N-s-b
      Y=dY[(s+1+b):N]
      X=dX[(s+1):(N-b)]
      if(s > 0){
         for (j in 1:s){
            X=cbind(X,dX[(s+1-j):(N-b-j)])
         }
      }
      m1=arima(Y,order=c(p,0,q),xreg=X)
      est=m1$coef; sigma2=m1$sigma2; residuals=m1$residuals; varcoef=m1$var.coef
      nx=dim(X)[2]
      se=sqrt(diag(m1$var.coef))
      coef.arma=NULL
      se.arma=NULL
      pq=p+q
      if(pq > 0){
         coef.arma=est[1:pq]
         se.arma=se[1:pq]
         p1=cbind(coef.arma,se.arma)
         cat("ARMA coefficients & s.e.:","\n")
         print(t(p1),digits=3)
      }
      v=est[(pq+1):(pq+1+nx)]
      se.v=se[(pq+1):(pq+1+nx)]
      pr=cbind(v,se.v)
      cat("Transfer function coefficients & s.e.:","\n")
      print(t(pr),digits=3)
   }
   else{
      ist=max(r,s)+1+b
      par=c(mean(dY[ist:N]))
      par=c(par,rep(0.1,s+1))
      par=c(par,rep(0.1,r))
      if(p > 0)par=c(par,rep(0.1,p))
      if(q > 0)par=c(par,rep(0.01,q))
      m11=nlm(Nlike,par,hessian=TRUE,dY=dY,dX=dX,orderN=orderN,orderX=orderX)
      est=m11$estimate
      varcoef=solve(m11$hessian)
      se=sqrt(diag(varcoef))
      residuals=Gaulike(est,dY=dY,dX=dX,orderN=orderN,orderX=orderX)
      sigma2=var(residuals)
      pq=p+q
      npar=length(est)
      v=est[1:(npar-pq)]
      se.v=se[1:(npar-pq)]
      pr=cbind(v,se.v)
      cat("Delay: ",b,"\n")
      cat("Transfer function coefficients & s.e.:","\n")
      cat("in the order: constant, omega, and delta:",c(1,s+1,r),"\n")
      print(t(pr),digits=3)
      if(pq > 0){
         coef.arma=est[(npar-pq+1):npar]
         se.arma=se[(npar-pq+1):npar]
         p1=cbind(coef.arma,se.arma)
         cat("ARMA order:","\n")
         print(c(p,dify,q))
         cat("ARMA coefficients & s.e.:","\n")
         print(t(p1),digits=3)
      }
      #
   }
   Nt = Nts(est,dY=dY,dX=dX,orderN=orderN,orderX=orderX)
   
   tfm1 <- list(estimate=est,sigma2=sigma2,residuals=residuals,varcoef=varcoef, Nt=Nt)
}
"tfm2" <- function(y,x,x2=NULL,ct=NULL,wt=NULL,orderN=c(1,0,0),orderS=c(0,0,0),sea=12,order1=c(0,1,0),order2=c(0,-1,0)){
   ## Estimation of a transfer function model with TWO exogenous variables
   ### The model is Y_t- c0 -c1*c_t -c2*w_t - w(B)/d(B)X_t  - W(b)/D(B)X_{2t} = theta(B)*Theta(B)/[phi(B)*Phi(B)]a_t.
   ### orderN = c(p,d,q) for the regular ARMA part
   ### orderS = c(P,D,Q) for the seasonal ARMA part
   ### order1 = c(r,s,b) where d(B) = 1 - d_1B - ... - d_r B^r
   ###            and w(B) = w_0+w_1B + ... + w_s B^s and b is the delay.
   ###
   ### order2 = c(r2,s2,b2) for the second exogenous variable
   ### wt: for co-integrated system 
   ### ct: a given determinsitic variable such as time trend
   ### par=c(c0,w0,w1,...,ws,d1,...,dr,c1,c2,W0, ...,Ws,D1,...,Dr,phi,theta,Phi,Theta): November 2014
   ###
   dify = orderN[2]; dY=y; dX=x; dX2=x2; dW=wt; dC=ct
   phi=NULL; theta=NULL; Phi=NULL; Theta=NULL; omega=NULL; delta=NULL; Omega=NULL; Delta=NULL
   if(dify > 0){
      dY <- y[(dify+1):length(y)]-y[1:(length(y)-dify)]
      dX <- x[(dify+1):length(x)]-x[1:(length(x)-dify)]
      if(!is.null(x2)){
        dX2 <- x2[(dify+1):length(x2)]-x2[1:(length(x2)-dify)]
        }
      if(!is.null(wt)){
        dW <- wt[(dify+1):length(wt)] - wt[1:(length(wt)-dify)]
        }
      if(!is.null(ct)){
       dC <- ct[(dify+1):length(ct)] - ct[1:(length(ct)-dify)]
       }
    }
### seasonal difference, if any
   difys = orderS[2]; lags=difys*sea
   if(difys > 0){
      dY <- dY[(lags+1):length(dY)]-dY[1:(length(dY)-lags)]
      dX <- dX[(lags+1):length(dX)]-dX[1:(length(dX)-lags)]
      if(!is.null(x2)){
        dX2 <- dX2[(lags+1):length(dX2)]-dX2[1:(length(dX2)-lags)]
        }
      if(!is.null(wt)){
        dW <- dW[(lags+1):length(dW)] - dW[1:(length(dW)-lags)]
        }
      if(!is.null(ct)){
       dC <- dC[(lags+1):length(dC)] - dC[1:(length(dC)-lags)]
        }
   }
#
   N = length(dY); N1=length(dX)
   N=min(N,N1)
   if(length(dX2) > 0)N=min(N,length(dX2))
   if(length(dW) > 0) N=min(N,length(dW))
   if(length(dC) > 0) N=min(N,length(dC))
   phi=NULL; theta=NULL; ome=NULL; del=NULL; ome2=NULL; del2=NULL; Phi=NULL; Theta=NULL
   r=order1[1]; s=order1[2]; b=order1[3]; p=orderN[1]; q=orderN[3]; P=orderS[1]; Q=orderS[3]
   r=max(r,0); s=max(0,s); b=max(0,b); p=max(0,p); q=max(0,q); P=max(0,P); Q=max(0,Q)
   r2=order2[1]; s2=order2[2]; b2=order2[3]
### subroutines used
 Nlike <- function(par,dY=dY,dX=dX,dX2=dX2,dW=dW,dC=dC,orderN=orderN,orderS=orderS,sea=sea,order1=order1,order2=order2){
   resi = Gaulike(par,dY=dY,dX=dX,dX2=dX2,dW=dW,dC=dC,orderN=orderN,orderS=orderS,sea=sea,order1=order1,order2=order2)
   sig=sqrt(var(resi))
   n1=length(resi)
   Nlike=-sum(dnorm(resi,mean=rep(0,n1),sd=sig,log=TRUE))
  }

 Gaulike <- function(par,dY=dY,dX=dX,dX2=dX2,dW=dW,dC=dC,orderN=orderN,orderS=orderS,sea=sea,order1=order1,order2=order2){
   p=orderN[1]; q=orderN[3]; r=order1[1]; s=order1[2]; b=order1[3]; P=orderS[1]; Q=orderS[3]
   r2=order2[2]; s2=order2[2]; b2=order2[3]
   c0=par[1]
   ome=par[2:(2+s)]
   #
   if(r > 0)del=par[(2+s+1):(2+s+r)]
   icnt=2+s+r
   if(!is.null(dC)){c1=par[icnt+1]
                   icnt=icnt+1
                   }
   if(!is.null(dW)){c2=par[icnt+1]
                   icnt=icnt+1
                   }
   if(!is.null(dX2)){
          ome2=par[(icnt+1):(icnt+1+s2)]
          icnt=icnt+1+s2
          if(r2 > 0){del2=par[(icnt+1):(icnt+r2)]
                      icnt=icnt+r2
                    }
          }
   if(p > 0){phi=par[(icnt+1):(icnt+p)]
             icnt=icnt+p
             }
   if(q > 0){theta=par[(icnt+1):(icnt+q)]
             icnt=icnt+q
            }
   if(P > 0){Phi=par[(icnt+1):(icnt+P)]
             icnt=icnt+P
            }
   if(Q > 0)Theta=par[(icnt+1):(icnt+Q)]
#
   N=length(dY)
   N1t=dY-c0
   Nt = dX
   if(r > 0){
     Nt=filter(dX,del,method="r",init=rep(mean(dX),r))
     }
   ##
   ist=max(b+s+1,b2+s2+1)
   N=length(Nt)
   N1t=N1t[ist:N]-ome[1]*Nt[(ist-b):(N-b)]
   if(s > 0){
      for (j in 1:s){
         N1t=N1t-ome[j+1]*Nt[(ist-j-b):(N-j-b)]
       }
    }
   if(!is.null(dC))N1t=N1t-c1*dC[ist:N]
   if(!is.null(dW))N1t=N1t-c2*dW[ist:N]
   if(!is.null(dX2)){
    Zt=dX2
    if(r2 > 0){
     Zt=filter(dX2,del2,method="r",init=rep(mean(dX2),r2))
     }
    N1t=N1t - ome2[1]*Zt[(ist-b2):(N-b2)]
    if(s2 > 0){
       for (j in 1:s2){
        N1t=N1t-ome2[j+1]*Zt[(ist-j-b2):(N-j-b2)]
        }
       }
    }
   N1=length(N1t)
   re=N1t[(p+1):N1]
   if(p > 0){
      for (j in 1:p){
         re=re-phi[j]*N1t[(p+1-j):(N1-j)]
        }
     }
   #
   if(q > 0)re=filter(re,theta,method="r",init=rep(0,q))
   N1=length(re)
   resi=re[(P*sea+1):N1]
   if(P > 0){
     for (j in 1:P){
      resi=resi-Phi[j]*re[(P*sea+1-j*sea):(N1-j*sea)]
      }
     }
    if(Q > 0){
     f1=rep(0,sea*Q)
     for (j in 1:Q){
      f1[j*sea]=Theta[j]
      }
     resi=filter(resi,f1,method="r",init=rep(0,sea*Q))
    }
   Gaulike = resi
 }

### Obtain the N(t) series
 Nts <- function(par,dY=dY,dX=dX,dX2=dX2,dW=dW,dC=dC,order1=order1,order2=order2){
   r=order1[1]; s=order1[2]; b=order1[3]
   r2=order2[1]; s2=order2[2]; b2=order2[3]
   c0=par[1]
   ome=par[2:(2+s)]
   icnt=2+s
   #
   if(r > 0){del=par[(2+s+1):(2+s+r)]
             icnt=2+s+r
            }
   if(!is.null(dC)){c1=par[icnt+1]
                   icnt=icnt+1
                  }
   if(!is.null(dW)){c2=par[icnt+1]
                   icnt=icnt+1
                 }
   if(!is.null(dX2)){
              ome2=par[(icnt+1):(icnt+1+s2)]
              icnt=icnt+1+s2
        if(r2 > 0){
             del2=par[(icnt+1):(icnt+r2)]
             icnt=icnt+r2
             }
        }
   N=length(dY)
   N1t=dY-c0
   Nt=dX
   if(r > 0){
     Nt=filter(dX,del,method="r",init=rep(mean(dX),r))
    }
   ##
   ist=max(b+s+1,b2+s2+1)
   N=length(Nt)
   N1t=N1t[ist:N]-ome[1]*Nt[(ist-b):(N-b)]
   if(s > 0){
      for (j in 1:s){
         N1t=N1t-ome[j+1]*Nt[(ist-j-b):(N-j-b)]
        }
     }
   if(!is.null(dC))N1t=N1t-c1*dC[ist:N]
   if(!is.null(dW)){N1t=N1t-c2*dW[ist:N]}
   if(!is.null(dX2)){
     Zt=dX2
     if(r2 > 0){
       Zt=filter(Zt,del2,method="r",init=rep(mean(dX2),r2))
       }
      N1t=N1t - ome2[1]*Zt[(ist-b2):(N-b2)]
      if(s2 > 0){
        for (j in 1:s2){
         N1t=N1t-ome2[j+1]*Zt[(ist-b2-j):(N-b2-j)]
         }
        }
     }
   Nts=N1t
  }
####
   ## r = 0 && r2=0, the model can be fitted by the regular "arima" command.
   if((r==0) && (r2==0)){
      ist=max(s+b,s2+b2)+1
      nobe=N-ist+1
      Y=dY[ist:N]
      X=dX[(ist-b):(N-b)]
      if(s > 0){
         for (j in 1:s){
            X=cbind(X,dX[(ist-b-j):(N-b-j)])
         }
      }
      if(!is.null(dC)){X=cbind(X,dC[ist:N])}
      if(!is.null(dW)){X=cbind(X,dW[ist:N])}
      if(!is.null(dX2)){
         X=cbind(X,dX2[(ist-b2):(N-b2)])
         if(s2 > 0){
           for (j  in 1:s2){
             X=cbind(X,dX2[(ist-b2-j):(N-b2-j)])
             }
            }
        }
      X=as.matrix(X)
      if(min(P,Q) > 0){
       m1=arima(Y,order=c(p,0,q),seasonal=list(order=c(P,0,Q),period=sea),xreg=X)
       }
       else{
        m1=arima(Y,order=c(p,0,q),xreg=X)
      }
     est=m1$coef; sigma2=m1$sigma2; residuals=m1$residuals; varcoef=m1$var.coef
#### Changing the sign of the MA coefficients, if any.
     if(q > 0){
      for (j in 1:q){
       loc=p+j
       est[loc]=-est[loc]
       }
      }
     if(Q > 0){
       for (j in 1:Q){
        loc=p+q+P+j
        est[loc]=-est[loc]
        }
      }
#### re-ordering the estimate for computing Nt series
     jcnt=p+q+P+Q+1
     est1=c(est[jcnt:(jcnt+ncol(X))],est[1:(jcnt-1)])
###   
      nx=dim(X)[2]
      se=sqrt(diag(m1$var.coef))
      coef.arma=NULL
      se.arma=NULL
      pq=p+q
      if(pq > 0){
         coef.arma=est[1:pq]
         se.arma=se[1:pq]
         p1=cbind(coef.arma,se.arma)
         cat("Regular ARMA coefficients & s.e.:","\n")
         print(t(p1),digits=3)
         if(p > 0)phi=coef.arma[1:p]
         if(q > 0)theta=coef.arma[(p+1):pq]
      }
      PQ=P+Q
      if(PQ > 0){
       coef.sea=est[(pq+1):(pq+PQ)]
       se.sea=se[(pq+1):(pq+PQ)]
       psea=cbind(coef.sea,se.sea)
       cat("Seasonal ARMA coefficients & s.e.: ","\n")
       print(t(psea),digits=3)
       if(P > 0)Phi=coef.sea[1:P]
       if(Q > 0)Theta=coef.sea[(P+1):PQ]
      }
      icnt=pq+PQ
      v=est[(icnt+1):(icnt+1+nx)]
      se.v=se[(icnt+1):(icnt+1+nx)]
      pr=cbind(v,se.v)
      cat("Transfer function coefficients & s.e.:","\n")
      print(t(pr),digits=3)
      cat("Sigma-square & sigma: ",c(sigma2,sqrt(sigma2)),"\n")
      omega=v[1:(s+1)]
      kcnt=s+1
      if(!is.null(dC))kcnt=kcnt+1
      if(!is.null(dW))kcnt=kcnt+1
      Omega=v[(kcnt+1):nx]
   est=est1
   }
   else{
      ist=max(r,s)+1+b
      ist1=max(r2,s2)+1+b2
      ist=max(ist,ist1)
      par=c(mean(dY[ist:N]))
      par=c(par,rep(0.1,s+1))
      par=c(par,rep(0.1,r))
      if(!is.null(dC))par=c(par,.01)
      if(!is.null(dW))par=c(par,.1)
      if(!is.null(dX2)){
       par=c(par,rep(0.1,s2+1))
       par=c(par,rep(0.1,r2))
       }
      if(p > 0)par=c(par,rep(0.1,p))
      if(q > 0)par=c(par,rep(0.01,q))
      if(P > 0)par=c(par,rep(0.01,P))
      if(Q > 0)par=c(par,rep(0.01,Q))
      m11=nlm(Nlike,par,hessian=TRUE,dY=dY,dX=dX,dX2=dX2,dW=dW,dC=dC,orderN=orderN,orderS=orderS,sea=sea,order1=order1,order2=order2)
      est=m11$estimate
      varcoef=solve(m11$hessian)
      se=sqrt(diag(varcoef))
      residuals=Gaulike(est,dY=dY,dX=dX,dX2=dX2,dW=dW,dC=dC,orderN=orderN,orderS=orderS,sea=sea,order1=order1,order2=order2)
      sigma2=var(residuals)
      pq=p+q
      PQ=P+Q
      icnt=1+s+1+r
      v=est[1:icnt]
      se.v=se[1:icnt]
      pr=cbind(v,se.v)
      cat("First exogenous variable: ","\n")
      cat("Delay: ",b,"\n")
      cat("Transfer function coefficients & s.e.:","\n")
      cat("in the order: constant, omega, and delta:",c(1,s+1,r),"\n")
      print(t(pr),digits=3)
      cnst=v[1]
      omega=v[2:(s+2)]
      if(r > 0)delta=v[(s+3):icnt]
      if(!is.null(dC)){icnt=icnt+1
       cat("co-integrated coefficient & se: ",c(est[icnt],se[icnt]),"\n")
       }
      if(!is.null(dW)){icnt=icnt+1
        cat("Co-integration coefficient & se: ",c(est[icnt],se[icnt]),"\n")
         }
      if(!is.null(dX2)){
       jcnt=1+s2+r2
       v=est[(icnt+1):(icnt+jcnt)]
       se.v=se[(icnt+1):(icnt+jcnt)]
       pr=cbind(v,se.v)
       cat("Second exogenous variable: ","\n")
       cat("Delay: ",b2,"\n")
       cat("The transfer function coefficients & s.e.:","\n")
       cat("in the order: omega2 and delta2: ",c(s2+1,r2),"\n")
       print(t(pr),digits=3)
       Omega=v[1:(s2+1)]
       if(r2 > 0)Delta=v[(2+s2):jcnt]
       icnt=icnt+jcnt
      }
      if(pq > 0){
         coef.arma=est[(icnt+1):(icnt+pq)]
         se.arma=se[(icnt+1):(icnt+pq)]
         p1=cbind(coef.arma,se.arma)
         cat("Regular ARMA order:","\n")
         print(c(p,dify,q))
         cat("Regular ARMA coefficients & s.e.:","\n")
         print(t(p1),digits=3)
         icnt=icnt+pq
         if(p > 0)phi=coef.arma[1:p]
         if(q > 0)theta=coef.arma[(p+1):pq]
        }
      if(PQ > 0){
        coef.sea=est[(icnt+1):(icnt+PQ)]
        se.sea=se[(icnt+1):(icnt+PQ)]
        ps=cbind(coef.sea,se.sea)
        cat("Seasonal ARMA order: ","\n")
        print(c(P,difys,Q))
        cat("Seasonal ARMA coefficients & s.e.: ","\n")
        print(t(ps),digits=3)
        if(P > 0)Phi=coef.sea[1:P]
        if(Q > 0)Theta=coef.sea[(P+1):PQ]
       }
    cat("Sigma-square & sigma: ",c(sigma2,sqrt(sigma2)),"\n")
   }
#
   Nt <- Nts(est,dY=dY,dX=dX,dX2=dX2,dW=dW,dC=dC,order1=order1,order2=order2)
   
   tfm2 <- list(estimate=est,sigma2=sigma2,residuals=residuals,varcoef=varcoef,Nt=Nt,rAR=phi,rMA=theta,sAR=Phi,sMA=Theta,
omega=omega,delta=delta,omega2=Omega,delta2=Delta)
}

### Back-testing
"Btfm2" <- function(y,x,x2=NULL,wt=NULL,ct=NULL,orderN=c(1,0,0),orderS=c(0,0,0),sea=12,order1=c(0,1,0),order2=c(0,-1,0),orig=(length(y)-1)){
err=NULL
r=order1[1]; s=order1[2]; b=order1[3]
r2=order2[1]; s2=order2[2]; b2=order2[3]
p = orderN[1]; dify=orderN[2]; q=orderN[3]
P = orderS[1]; difys=orderS[2]; Q=orderS[3]
dY=y; dX=x; dX2=x2; dW=wt; dC=ct
if(dify > 0){
      dY <- y[(dify+1):length(y)]-y[1:(length(y)-dify)]
      dX <- x[(dify+1):length(x)]-x[1:(length(x)-dify)]
      if(!is.null(x2)){
        dX2 <- x2[(dify+1):length(x2)]-x2[1:(length(x2)-dify)]
        }
      if(!is.null(wt)){
        dW <- wt[(dify+1):length(wt)] - wt[1:(length(wt)-dify)]
        }
      if(!is.null(ct)){
        dC <- ct[(dify+1):length(ct)] - ct[1:(length(ct)-dify)]
       }
    }
if(difys > 0){
      lags=difys*sea
      dY <- dY[(lags+1):length(dY)]-dY[1:(length(dY)-lags)]
      dX <- dX[(lags+1):length(x)]-dX[1:(length(dX)-lags)]
      if(!is.null(x2)){
        dX2 <- dX2[(lags+1):length(dX2)]-dX2[1:(length(dX2)-lags)]
        }
      if(!is.null(wt)){
        dW <- dW[(lags+1):length(dW)] - dW[1:(length(dW)-lags)]
        }
      if(!is.null(ct)){
        dC <- dC[(lags+1):length(dC)] - dC[1:(length(dC)-lags)]
       }
    }
   N = length(dY); N1=length(dX)
   N=min(N,N1)
   if(length(dX2) > 0)N=min(N,length(dX2))
   if(length(dW) > 0) N=min(N,length(dW))
   if(length(dC) > 0) N=min(N,length(dC))
   orig=orig-dify-difys*sea
## function to perform prediction: 1-step ahead only
###
fore1 <- function(par,dY=dY,dX=dX,dX2=x2p,dW=wtp,dC=ctp,orderN=orderN,orderS=orderS,sea=sea,order1=order1,order2=order2,resi=resi){
   p=orderN[1]; q=orderN[3]; r=order1[1]; s=order1[2]; b=order1[3]
   r2=order2[2]; s2=order2[2]; b2=order2[3]; P=orderS[1]; Q=orderS[3]
   c0=par[1]
   ome=par[2:(2+s)]
   #
   if(r > 0)del=par[(2+s+1):(2+s+r)]
   icnt=2+s+r
   if(!is.null(dC)){c1=par[icnt+1]
                  icnt=icnt+1
                  }
   if(!is.null(dW)){c2=par[icnt+1]
                   icnt=icnt+1
                   }
   if(!is.null(dX2)){
          ome2=par[(icnt+1):(icnt+1+s2)]
          icnt=icnt+1+s2
          if(r2 > 0){del2=par[(icnt+1):(icnt+r2)]
                      icnt=icnt+r2
                    }
          }
   if(p > 0){phi=par[(icnt+1):(icnt+p)]
             icnt=icnt+p
             }
   if(q > 0){theta=par[(icnt+1):(icnt+q)]
            icnt=icnt+q
            }
   if(P > 0){Phi=par[(icnt+1):(icnt+P)]
             icnt=icnt+P
            }
   if(Q > 0){Theta=par[(icnt+1):(icnt+Q)]
            }
   N=length(dY)
   tmp=dY-c0
   Nt = dX
   if(r > 0){
     Nt=filter(dX,del,method="r",init=rep(mean(dX),r))
     }
   if(b == 0){
     tmp=tmp-ome[1]*Nt
     }else{
     tmp=tmp-ome[1]*c(rep(0,b),Nt[1:(N-b)])
     }
   if(s > 0){
    for (j in 1:s){
     tmp=tmp-ome[j+1]*c(rep(0,b+j),Nt[1:(N-b-j)])
     }
    }
   if(!is.null(dC))tmp=tmp-c1*dC
   if(!is.null(dW))tmp=tmp-c2*dW
   if(!is.null(dX2)){
    Zt=dX2
    if(r2 > 0){
     Zt=filter(dX2,del2,method="r",init=rep(mean(dX2),r2))
     }
    tmp=tmp-ome2[1]*c(rep(0,b2),Zt[1:(N-b2)])
    if(s2 > 0){
       for (j in 1:s2){
        tmp=tmp-ome2[j+1]*c(rep(0,j+b2),Zt[1:(N-j-b2)])
        }
       }
    }
### The next step is for prediction, starting with exogenous variables at time t+1.
   pred=dY[N]-tmp[N]
   if(p > 0){
    for (j in 1:p){
     pred=pred+phi[j]*tmp[N-j]
     }
    }
   if(P > 0){
     for (j in 1:P){
      pred=pred+Phi[j]*tmp[N-j*sea]
      }
     }
   if((p > 0)&&(P > 0)){
     for (j in 1:P){
       j1=j*sea
       for (i in 1:p){
        pred=pred-Phi[j]*phi[i]*tmp[N-j1-i]
        }
       }
     }
   if(q > 0){
     for (j in 1:q){
      pred=pred-theta[j]*resi[length(resi)+1-j]
      }
     }
    if(Q > 0){
     for (j in 1:Q){
       pred=pred-Theta[j]*resi[length(resi)+1-j*sea]
       }
      }
   if((q > 0)&&(Q > 0)){
     for (j in 1:Q){
      j1=j*sea
      for (i in 1:q){
       pred=pred+Theta[j]*theta[i]*resi[length(resi)+1-i-j1]
       }
      }
    }
   err=dY[N]-pred
   err
}
###
nT=length(dY)
if(nT > orig){
### Estimation
for (it in orig:(nT-1)){
  x2p=NULL; wtp=NULL; ctp = NULL
  if(!is.null(x2))x2p=dX2[1:it]
  if(!is.null(wt))wtp=dW[1:it]
  if(!is.null(ct))ctp=dC[1:it]
  m1 = tfm2(dY[1:it],dX[1:it],x2=x2p,wt=wtp,ct=ctp,orderN=orderN,orderS=orderS,sea=sea,order1=order1,order2=order2)
  par=m1$estimate
  Tp1=it+1
  resi=m1$residuals
  nr=length(resi)
  if(nr < it){
   resi=c(rep(0,it-nr),resi)
  }
### prediction via computing the residuals
 x2p=NULL; wtp=NULL; ctp=NULL
 if(!is.null(x2))x2p=dX2[1:Tp1]
 if(!is.null(wt))wtp=dW[1:Tp1]
 if(!is.null(ct))ctp=dC[1:Tp1]
 error = fore1(par,dY=dY[1:Tp1],dX=dX[1:Tp1],dX2=x2p,dW=wtp,dC=ctp,orderN=orderN,orderS=orderS,sea=sea,order1=order1,order2=order2,resi=resi)
##
 err=c(err,error)
 }
}
 bias=mean(err); nf=length(err)
 mse=mean(err^2); mae=mean(abs(err))
 rmse=sqrt(mse)
 cat("Forecast origin & number of forecasts: ",c(orig,nf),"\n")
 cat("bias,  mse, rmse & MAE: ",c(bias, mse,rmse, mae),"\n")
 Btfm2 <- list(ferror=err,mse=mse,rmse=rmse,mae=mae,nobf=nf)
}

####
"VARchi" <- function(x,p=1,include.mean=T,thres=1.645){
   # Fits a vector AR(p) model, then performs
   # a chi-square test to zero out insignificant parameters.
   if(!is.matrix(x))x=as.matrix(x)
   Tn=dim(x)[1]
   k=dim(x)[2]
   if(p < 1)p=1
   ne=Tn-p
   ist=p+1
   y=x[ist:Tn,]
   if(include.mean){
      xmtx=cbind(rep(1,ne),x[p:(Tn-1),])
   }
   else {
      xmtx=x[p:(Tn-1),]
   }
   if(p > 1){
      for (i in 2:p){
         xmtx=cbind(xmtx,x[(ist-i):(Tn-i),])
      }
   }
   #
   #perform estimation
   ndim=dim(xmtx)[2]
   res=NULL
   xm=as.matrix(xmtx)
   xpx=crossprod(xm,xm)
   xpxinv=solve(xpx)
   xpy=t(xm)%*%as.matrix(y)
   beta=xpxinv%*%xpy
   resi=y-xm%*%beta
   sse=t(resi)%*%resi/(Tn-p-ndim)
   C1=kronecker(sse,xpxinv)
   dd=sqrt(diag(C1))
   #
   bhat=c(beta)
   tratio=bhat/dd
   para=cbind(bhat,dd,tratio)
   npar=length(bhat)
   K=NULL
   omega=NULL
   for (i in 1:npar){
      if(abs(tratio[i]) < thres){
         idx=rep(0,npar)
         idx[i]=1
         K=rbind(K,idx)
         omega=c(omega,bhat[i])
      }
   }
   v=dim(K)[1]
   K=as.matrix(K)
   cat("Number of targeted parameters: ",v,"\n")
   #####print(K)
   if(v > 0){
      C2=K%*%C1%*%t(K)
      C2inv=solve(C2)
      tmp=C2inv%*%as.matrix(omega,v,1)
      chi=sum(omega*tmp)
      pvalue=1-pchisq(chi,v)
      cat("Chi-square test and p-value: ",c(chi,pvalue),"\n")
   }
   else{
      print("No contraints needed")
   }
   VARchi<-list(data=x,cnst=include.mean,order=p,coef=beta,constraints=K,omega=omega,covomega=C2)
}

###
"FEVdec" <- function(Phi,Theta,Sig,lag=4){
   # Perform forecast error vcovariance decomposition
   #
   # Phi: k by kp matrix of AR coefficients, i.e. [AR1,AR2,AR3, ..., ARp]
   # Theta: k by kq matrix of MA coefficients, i.e. [MA1,MA2, ..., MAq]
   # Sig: residual covariance matrix
   # Output: (a) Plot and (b) Decomposition
   if(length(Phi) > 0){
      if(!is.matrix(Phi))Phi=as.matrix(Phi)
     }
   if(length(Theta) > 0){
      if(!is.matrix(Theta))Theta=as.matrix(Theta)
     }
   if(!is.matrix(Sig))Sig=as.matrix(Sig)
   if(lag < 1) lag=1
   # Compute MA representions: This gives impulse response function without considering Sigma.
   p = 0
   if(length(Phi) > 0){
      k=nrow(Phi)
      m=ncol(Phi)
      p=floor(m/k)
    }
   q=0
   if(length(Theta) > 0){
      k=dim(Theta)[1]
      m=dim(Theta)[2]
      q=floor(m/k)
    }
   cat("Order of the ARMA mdoel: ","\n")
   print(c(p,q))
   # Consider the MA part to psi-weights
   Si=diag(rep(1,k))
   if(q > 0){
      Si=cbind(Si,-Theta)
    }
   m=(lag+1)*k
   m1=(q+1)*k
   if(m > m1){
      Si=cbind(Si,matrix(0,k,(m-m1)))
    }
   #
   if (p > 0){
      for (i in 1:lag){
         if (i <= p){
            idx=(i-1)*k
            tmp=Phi[,(idx+1):(idx+k)]
         }
         else{
            tmp=matrix(0,k,k)
         }
         #
         jj=i-1
         jp=min(jj,p)
         if(jp > 0){
            for(j in 1:jp){
               jdx=(j-1)*k
               idx=(i-j)*k
               w1=Phi[,(jdx+1):(jdx+k)]
               w2=Si[,(idx+1):(idx+k)]
               tmp=tmp+w1%*%w2
               ##print(tmp,digits=4)
            }
         }
         kdx=i*k
         Si[,(kdx+1):(kdx+k)]=tmp
         ## end of i loop
      }
      ## end of (p > 0)
     }
   # Compute the impulse response of orthogonal innovations
   orSi=NULL
   m1=chol(Sig)
   P=t(m1)
   orSi=P
   for(i in 1:lag){
      idx=i*k
      w1=Si[,(idx+1):(idx+k)]
      w2=w1%*%P
      orSi=cbind(orSi,w2)
   }
   #### Compute the covariance matrix of forecast errors
   orSi2=orSi^2
   ##### compute the partial sum (summing over lags)
   Ome=orSi2[,1:k]
   wk=Ome
   for (i in 1:lag){
      idx=i*k
      wk=wk+orSi2[,(idx+1):(idx+k)]
      Ome=cbind(Ome,wk)
   }
   FeV=NULL
   ##
   OmeRa = Ome[,1:k]
   FeV=cbind(FeV,apply(OmeRa,1,sum))
   OmeRa = OmeRa/FeV[,1]
   for (i in 1:lag){
      idx=i*k
      wk=Ome[,(idx+1):(idx+k)]
      FeV=cbind(FeV,apply(wk,1,sum))
      OmeRa=cbind(OmeRa,wk/FeV[,(i+1)])
     }
   cat("Standard deviation of forecast error: ","\n")
   print(sqrt(FeV))
   #
   cat("Forecast-Error-Variance Decomposition","\n")
   for (i in 1:(lag+1)){
      idx=(i-1)*k
      cat("Forecast horizon: ",i,"\n")
      Ratio=OmeRa[,(idx+1):(idx+k)]
      print(Ratio)
    }
   FEVdec <- list(irf=Si,orthirf=orSi,Omega=Ome,OmegaR=OmeRa)
  }

####
"mFilter" <- function(da,Wgt,init=NULL){
   # Multivariate filtering algorithm: using the nagative pi-weights
   ## (i+pi1 B + pi2 B^2 + ....)z_t = a_t.
   # Created by Ruey S. Tsay in April 2012.
   #
   # Wgt=[Theta1, Theta2,..., Thetaq]
   # Filered data = a_t - Theta1 *a_{t-1} - ... - Thetaq * a_{t-q}.
   #
   if(!is.matrix(da))da=as.matrix(da)
   if(!is.matrix(Wgt))Wgt=as.matrix(Wgt)
   if(length(init) > 0)init=as.matrix(init)
   #
   #### set up the data matrix for filtering
   ### In mFilter
   ##cat("in mFilter: ","\n")
   ##print(Wgt)
   
   nT=dim(da)[1]
   k=dim(da)[2]
   if(k == 1){
      q=length(Wgt)
      Wgt=matrix(Wgt,1,q)
   }
   else{
      m=dim(Wgt)[2]
      q=floor(m/k)
   }
   if(length(init) < 0){
      nit = 0
      x=da}
   else{
      nit=dim(init)[1]
      x=rbind(init,da)
   }
   if(k == 1) x=matrix(x,length(x),1)
   # obtain the nagative pi-weights
   Npi = diag(rep(1,k))
   Npi = cbind(Wgt[,1:k],Npi)
   TT=dim(x)[1]
   ####
   for (i in 2:(TT-1)){
      kend=min(q,i)
      if(k == 1){
         tmp=0
         for (j in 1:kend){
            jdx=j-1
            w1=Npi[jdx+1]
            w2=Wgt[jdx+1]
            tmp=tmp+w1*w2
         }
         Npi=c(tmp,Npi)
      }
      else{
         tmp=matrix(0,k,k)
         for (j in 1:kend){
            jdx=(j-1)*k
            w1=Npi[,(jdx+1):(jdx+k)]
            w2=Wgt[,(jdx+1):(jdx+k)]
            tmp=tmp+w1%*%w2
         }
         Npi=cbind(tmp,Npi)
      }
      #
   }
   T1=dim(Npi)[2]
   Wrk=c(t(x))
   T2=length(Wrk)
   At=NULL
   for (it in 1:nT){
      mk=(it-1)*k
      wk3=Npi[,(mk+1):T1]%*%Wrk[1:(T2-mk)]
      At=rbind(t(wk3),At)
      ##if(it < 4)print(At)
   }
   At
 }


#### Exact likelihood VMA programs
"VMAe" <- function(da,q=1,include.mean=T,coef0=NULL,secoef0=NULL,fixed=NULL,prelim=F,details=F,thres=2.0){
   # Estimation of a vector MA model using EXACT MLE (Gaussian dist)
   ### coef0 and secoef0 are the initial estimates and their standard errors (mainly from the conditional estimates).
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
   #xtx=t(xmtx)%*%xmtx
   xtx=crossprod(xmtx,xmtx)
   #xty=t(xmtx)%*%ymtx
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
##
   if(length(coef0) < 1){
      # Obtain initial parameter estimates if necessary
      ### Use VAR approximation to obtain initial parameter estimates
      m1=VARorder(da,q+12,output=FALSE)
      porder=m1$aicor
      if(porder < 1)porder=1
      m2=VAR(da,porder,output=FALSE)
      y=da[(porder+1):nT,]
      x=m2$residuals
      m3=THini(y,x,q,include.mean)
      beta=-m3$estimates
      sebeta=m3$se
      nr=dim(beta)[1]
      if(include.mean){
         beta[1,]=-beta[1,]
      }
      ### Preliminary simplification
      if(prelim){
         fixed = matrix(0,nr,k)
         for (j in 1:k){
            tt=beta[,j]/sebeta[,j]
            idx=c(1:nr)[abs(tt) >= thres]
            fixed[idx,j]=1
         }
      }
      ## end initial estimation
   }
   else {
      beta=coef0
      sebeta=secoef0
      nr=dim(beta)[1]
   }
   #
   if(length(fixed)==0){fixed=matrix(1,nr,k)}
   #
   par=NULL
   separ=NULL
   #
   VMAecnt = 0
   ist=0
   if(include.mean){
      jdx=c(1:k)[fixed[1,]==1]
      VMAecnt=length(jdx)
      if(VMAecnt > 0){
         par=beta[1,jdx]
         separ=sebeta[1,jdx]
      }
      TH=beta[2:(kq+1),]
      seTH=sebeta[2:(kq+1),]
      ist=1
   }
   else {
      TH=beta
      seTH=sebeta
   }
   #########
   for (j in 1:k){
      idx=c(1:(nr-ist))[fixed[(ist+1):nr,j]==1]
      if(length(idx) > 0){
         par=c(par,TH[idx,j])
         separ=c(separ,seTH[idx,j])
      }
   }
   ###
   ParE <- par
   cat("Number of parameters: ",length(par),"\n")
   cat("initial estimates: ",par,"\n")
   ### Set up lower and upper bounds
   lowerBounds=par; upperBounds=par
   npar=length(par)
   mult=2.0
   if((npar > 10)||(q > 2))mult=1.5
   if(length(coef0) > 0){
      mult=1.0
   }
   #
   for (j in 1:npar){
      lowerBounds[j] = par[j]-mult*separ[j]
      upperBounds[j] = par[j]+mult*separ[j]
   }
   cat("Par. Lower-bounds: ",lowerBounds,"\n")
   cat("Par. Upper-bounds: ",upperBounds,"\n")
### likelihood function
 EVMAq <- function(par,zt=da,q=q,include.mean=include.mean,fixed=fixed,EstStep=T){
   # The model used is
   ## a_t' = x_t' - mu' + a_{t-1}'theta_1'+a_{t-2}'theta_2' + ....
   k=dim(zt)[2]
   nT=dim(zt)[1]
   #
   mu=rep(0,k)
   icnt=0; VMAecnt <- 0
   ist=0
   if(include.mean){
      ist=1
      jdx=c(1:k)[fixed[1,]==1]
      icnt=length(jdx); VMAecnt <- icnt
      if(icnt > 0)
       mu[jdx]=par[1:icnt]
   }
   ### remove the mean
   for (j in 1:k){
      zt[,j]=zt[,j]-mu[j]
   }
   ## obtain the Theta-matrix
   kq=k*q
   theta=matrix(0,kq,k)
   for (j in 1:k){
      idx=c(1:kq)[fixed[(ist+1):(ist+kq),j]==1]
      jcnt=length(idx)
      if(jcnt > 0){
         theta[idx,j]=par[(icnt+1):(icnt+jcnt)]
         icnt=icnt+jcnt
      }
   }
   # theta = rbind[theta_1',theta_2', ..., theta_q']
   theta=t(theta)
   ### Check for invertibility before applying mFilter
   ### If necessary, set up the expanded VMA(1) model and data.
   k1=dim(theta)[2]
   Theta=theta[,1:k]
   Zt=zt
   est=NULL
   if(q > 1){
      z0=cbind(diag(1,k*(q-1)),matrix(0,k*(q-1),k))
      Theta=rbind(theta,z0)
      Zt=cbind(zt,matrix(0,nT,k*(q-1)))
   }
   m1=eigen(Theta)
   V1=m1$values
   M1=Mod(V1)
   ich=0
   for (i in 1:k1){
      if(M1[i] > 1){
         V1[i]=1/V1[i]
         ich=1
      }}
   if(ich > 0){
      ###cat("Eigenvalue detection occurred","\n")
      P1=m1$vectors
      P1i=solve(P1)
      Theta=Re(P1%*%diag(V1)%*%P1i)
      ##print(Theta)
      ## Replace the re-normalized parameters
      beta=t(Theta[1:k,])
      ist=0; nr=kq;
      if(include.mean){
         ist=1
         nr=kq+1
      }
      #########
      for (j in 1:k){
         idx=c(1:k1)[fixed[(ist+1):nr,j]==1]
         if(length(idx) > 0){
            est=c(est,beta[idx,j])
         }
      }
 ##     par=est
   }
   #### DO not need to repace parameter estimate in evaluating Hessian
   if(EstStep){
      if(VMAecnt > 0){
         par = c(par[1:VMAecnt],est)
       }
      else {
         par = est
       }
    }
   ##
   theta=Theta[1:k,]
   at=mFilter(zt,theta)
   sig=t(at)%*%at/nT
   ## Obtain the square-root matrix of Sigma
   sig=(sig+t(sig))/2
   m1=eigen(sig)
   va=m1$values+10^(-10)
   VA=diag(1/sqrt(va))
   P=m1$vectors
   SigH=P%*%VA%*%t(P)
   if(q > 1)SigH=kronecker(diag(rep(1,q)),SigH)
   ##### k1 is the dimension of the expanded VMA(1) model when q > 1.
   k1= k*q
   #### Obtain the pi-wights (negative) for VMA(1) model
   ##### In the process, also obtain the X-tilde matrix
   Psi=diag(rep(1,k1))
   Psi=cbind(Theta,Psi)
   tmp=Theta
   X=-SigH
   tmp1=-SigH%*%tmp
   X=rbind(X,tmp1)
   if(nT > 2){
      for (i in 2:nT){
         tmp=tmp%*%Theta
         Psi=cbind(tmp,Psi)
         tmp1=-SigH%*%tmp
         X=rbind(X,tmp1)
      }
      # end of the statement if(nT > 2)
   }
   ## Obtain the intial estimate
   vZt=c(t(Zt))
   Y=rep(0,k1)
   nPsi=dim(Psi)[2]
   for (it in 1:nT){
      iend=it*k1
      wk1=vZt[1:iend]
      wk2=Psi[,(nPsi-iend+1):nPsi]
      wk=wk2%*%as.matrix(wk1,iend,1)
      Y=c(Y,SigH%*%wk)
   }
   XpX=crossprod(X,X)
   XpY=crossprod(X,Y)
   a0H=solve(XpX,XpY)
   resi=Y-X%*%a0H
   SSr=sum(resi^2)
   d1=det(XpX)
   d2=det(sig)
   llike=0.5*(SSr + nT*log(2*pi*d2) + log(d2))
   llike
  }

  # Step 5: Estimate Parameters and Compute Numerically Hessian:
   if(details){
      fit = nlminb(start = ParE, objective = EVMAq, zt=da,q=q,include.mean=include.mean,fixed=fixed, 
                lower = lowerBounds, upper = upperBounds, control = list(trace=3))}
   else{
      fit = nlminb(start = ParE, objective = EVMAq, zt=da,q=q,include.mean=include.mean,fixed=fixed,
       control=list(step.min=0.4,step.max=0.8), lower = lowerBounds, upper = upperBounds)
     }
   #
   est=fit$par
   ###
   ################### Checking for invertibility of the fitted VMA models.
   zt = da
   ist=0
   mu=rep(0,k)
   icnt=0
   if(include.mean){
      ist=1
      jdx=c(1:k)[fixed[1,]==1]
      icnt=length(jdx)
      if(icnt > 0)
      mu[jdx]=est[1:icnt]
   }
   ### remove the mean
   for (j in 1:k){
      zt[,j]=zt[,j]-mu[j]
   }
   ## obtain the Theta-matrix
   theta=matrix(0,kq,k)
   for (j in 1:k){
      idx=c(1:kq)[fixed[(ist+1):(ist+kq),j]==1]
      jcnt=length(idx)
      if(jcnt > 0){
         theta[idx,j]=est[(icnt+1):(icnt+jcnt)]
         icnt=icnt+jcnt
      }
   }
   # theta = rbind[theta_1',theta_2', ..., theta_q']
   theta=t(theta)
   ### Check for invertibility of the final estimates
   ### If necessary, set up the expanded VMA(1) model and data.
   k1=dim(theta)[2]
   Theta=theta[,1:k]
   Zt=zt
   if(q > 1){
      z0=cbind(diag(1,k*(q-1)),matrix(0,k*(q-1),k))
      Theta=rbind(theta,z0)
      Zt=cbind(zt,matrix(0,nT,k*(q-1)))
   }
   m1=eigen(Theta)
   V1=m1$values
   M1=Mod(V1)
   ich=0
   for (i in 1:k1){
      if(M1[i] > 1){
         V1[i]=1/V1[i]
         ich=1
      }}
   est1=NULL
   if(VMAecnt >0) est1=est[1:VMAecnt]
   if(ich > 0){
      ###cat("Eigenvalue detection occurred","\n")
      P1=m1$vectors
      P1i=solve(P1)
      Theta=Re(P1%*%diag(V1)%*%P1i)
      ##print(Theta)
      ## Replace the re-normalized parameters
      beta=t(Theta[1:k,])
      ist=0; nr=kq
      if(include.mean){
         ist=1
         nr=kq+1
      }
      #########
      for (j in 1:k){
         idx=c(1:k1)[fixed[(ist+1):nr,j]==1]
         if(length(idx) > 0){
            est1=c(est1,beta[idx,j])
         }
      }
      est=est1
   }
   ### The above steps of checking invertibility were added on April 22, 2012.
   cat("Final   Estimates: ",est,"\n")
   epsilon = 0.0001 * est
   npar=length(par)
   Hessian = matrix(0, ncol = npar, nrow = npar)
   for (i in 1:npar) {
      for (j in 1:npar) {
         x1 = x2 = x3 = x4  = fit$par
         x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
         x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
         x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
         x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
         Hessian[i, j] = (EVMAq(x1,zt=da,q=q,include.mean=include.mean,fixed=fixed,EstStep=F)
                         -EVMAq(x2,zt=da,q=q,include.mean=include.mean,fixed=fixed,EstStep=F)
                         -EVMAq(x3,zt=da,q=q,include.mean=include.mean,fixed=fixed,EstStep=F)
                         +EVMAq(x4,zt=da,q=q,include.mean=include.mean,fixed=fixed,EstStep=F))/
         (4*epsilon[i]*epsilon[j])
      }
   }
   # Step 6: Create and Print Summary Report:
   se.coef = sqrt(diag(solve(Hessian)))
   tval = est/se.coef
   matcoef = cbind(est, se.coef, tval, 2*(1-pnorm(abs(tval))))
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
      jdx=c(1:k)[fixed[1,]==1]
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
      idx=c(1:kq)[fixed[(ist+1):nr,j]==1]
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
   Past=matrix(0,1,kq)
   at=NULL
   for (t in 1:nT){
      tmp=zt[t,]+Past%*%TH
      at=rbind(at,tmp)
      if(q==1){
         Past=tmp
      }
      else{
         Past=c(tmp,Past[1:(kq-k)])
      }
   }
   sig=t(at)%*%at/nT
   cat(" ","\n")
   cat("Residuals cov-matrix:","\n")
   print(sig)
   dd=det(sig)
   d1=log(dd)
   aic=d1+2*npar/nT
   bic=d1+log(T)*npar/nT
   cat("----","\n")
   cat("aic= ",aic,"\n")
   cat("bic= ",bic,"\n")
   ### prepare fot output storage
   Theta=t(TH)
   if(include.mean){
      TH=rbind(cnt,TH)
      seTH=rbind(secnt,seTH)
   }
   
   VMAe <- list(data=da,MAorder=q,cnst=include.mean,coef=TH,secoef=seTH,residuals=at,Sigma=sig,Theta=Theta,mu=cnt,aic=aic,bic=bic)
 }

"refVMAe" <- function(model,thres=1){
   # This program refines the fitted models of VMA output by removing
   # insigificant parameters with abs(t-ratio) < thres.
   # model: output object from VMA
   x = model$data
   q = model$MAorder
   cnst = model$cnst
   coef=as.matrix(model$coef)
   secoef=as.matrix(model$secoef)
   nr=dim(coef)[1]
   nc=dim(coef)[2]
   for (j in 1:nc){
      for (i in 1:nr){
         if(secoef[i,j] < 10^(-8))secoef[i,j]=1.0
      }
    }
   fix=matrix(0,nr,nc)
   for (j in 1:nc){
      tt=coef[,j]/secoef[,j]
      idx=c(1:nr)[abs(tt) >= thres]
      fix[idx,j]=1
     }
   if(cnst){
      tt=coef[1,]/secoef[1,]
      idx=c(1:nc)[abs(tt) > 1.0]
      if(length(idx) > 0)fix[1,idx]=1
      }
     mm=VMAe(x,q=q,include.mean=cnst,fixed=fix,coef0=coef,secoef0=secoef)
    refVMAe <- list(data=x,MAorder=q,cnst=cnst,coef=mm$coef,secoef=mm$secoef,residuals=mm$residuals,Sigma=mm$Sigma,aic=mm$aic,bic=mm$bic,mu=mm$mu,Theta=mm$Theta)   
  }

###
"diffM" <- function(zt,d=1){
   ## taking difference of a vector time series
   ## d: (1-B^d)
   if(!is.matrix(zt))zt=as.matrix(zt)
   nT=dim(zt)[1]
   dzt = zt[(d+1):nT,]-zt[1:(nT-d),]
   dzt
}

##### Psi-weight calculation for a VARMA model.
"PSIwgt" <- function(Phi=NULL,Theta=NULL,lag=12,plot=TRUE,output=FALSE){
   ### Compute the psi-weight matrices of a VARMA(p,q) model,
   #### Phi=[phi1, phi2, ..., phip]
   #### Theta=[theta1,theta2,...,thetaq]
   #### Sigma= residual covariance matrix
   q=0; p=0; k=0
   if(length(Theta) > 0){
      k=dim(Theta)[1]
      k1=dim(Theta)[2]
      q=floor(k1/k)
   }
   #
   if(length(Phi) > 0){
      k=dim(Phi)[1]
      k1=dim(Phi)[2]
      p=floor(k1/k)
   }
   #
   if(k < 1) k=1
   PSI=diag(k); WGT=c(PSI)
   #
   for (il in 1:lag){
      ilk=il*k
      tmp=matrix(0,k,k)
      if((q > 0) && (il <= q))tmp=-Theta[,(ilk-k+1):ilk]
      if(p > 0){
         iend=min(il,p)
         for (j in 1:iend){
            jdx=(il-j)
            kdx=j*k
            tmp=tmp+Phi[,(kdx-k+1):kdx]%*%PSI[,(jdx*k+1):(jdx*k+k)]
         }
         ## end  of p > 0.
      }
      PSI=cbind(PSI,tmp)
      WGT=cbind(WGT,c(tmp))
      ### end of il-loop
   }
   ## print the output if needed
   if(output){
      for (i in 1:lag){
         cat("Lag: ",i," psi-matrix","\n")
         ist=i*k
         print(round(PSI[,(ist+1):(ist+k)],5))
      }
      ## end print
   }
   ## plots the psi-weights
   if(plot){
      tdx=c(1:(lag+1))-1
      par(mfcol=c(k,k),mai=c(0.3,0.3,0.3,0.3))
      gmax=max(WGT)
      gmin=min(WGT)
      cx=(gmax-gmin)/10
      gmax=gmax+cx
      gmin=gmin-cx
      for(j in 1:k^2){
         plot(tdx,WGT[j,],type='l',xlab='lag',ylab='Psiwgt',ylim=c(gmin,gmax),cex.axis=0.8)
         points(tdx,WGT[j,],pch='*',cex=0.8)
         title(main="Psi-weights")
      }
     par(mfcol=c(1,1))
    }
    PSIwgt <- list(psi.weight=PSI,irf=WGT)
  }

##### Pi-weight calculation for a VARMA model.
"PIwgt" <- function(Phi=NULL,Theta=NULL,lag=12,plot=TRUE){
   ### Compute the psi-weight matrices of a VARMA(p,q) model,
   #### Phi=[phi1, phi2, ..., phip]
   #### Theta=[theta1,theta2,...,thetaq]
   #### Sigma= residual covariance matrix
   m1=PSIwgt(Phi=Theta,Theta=Phi,lag=lag,plot=FALSE)
   PImtx=m1$psi.weight
   ###print(PImtx)
   k=dim(PImtx)[1]; nc=dim(PImtx)[2]
   PImtx[,(k+1):nc]=-PImtx[,(k+1):nc]
   lag=floor(nc/k)-1
   WGT=c(diag(k))
   for (i in 1:lag){
      cat("Lag: ",i," pi-matrix","\n")
      ist=(i-1)*k
      WGT=cbind(WGT,c(PImtx[,(ist+1):(ist+k)]))
      print(round(PImtx[,(ist+1):(ist+k)],5))
   }
   ## plots the pi-weights
   if(plot){
      tdx=c(1:(lag+1))-1
      par(mfcol=c(k,k),mai=c(0.3,0.3,0.3,0.3))
      gmax=max(WGT)
      gmin=min(WGT)
      cx=(gmax-gmin)/10
      gmax=gmax+cx
      gmin=gmin-cx
      for(j in 1:k^2){
         plot(tdx,WGT[j,],type='l',xlab='lag',ylab='Piwgt',ylim=c(gmin,gmax),cex.axis=0.8)
         points(tdx,WGT[j,],pch='*',cex=0.8)
         title(main="Pi-weights")
       }
     par(mfcol=c(1,1))
     }
   PIwgt <- list(pi.weight=PImtx)
 }

"VARMAirf" <- function(Phi=NULL,Theta=NULL,Sigma=NULL,lag=12,orth=TRUE){
   #### Phi=[phi1, phi2, ..., phip]
   #### Theta=[theta1,theta2,...,thetaq]
    q=0; p=0; k=0
   if(length(Theta) > 0){
      k=dim(Theta)[1]
      k1=dim(Theta)[2]
      q=floor(k1/k)
   }
   #
   if(length(Phi) > 0){
      k=dim(Phi)[1]
      k1=dim(Phi)[2]
      p=floor(k1/k)
   }
   #
   if(is.null(Sigma)){
      Sigma=diag(rep(1,k))
   }
   #
   if(orth){
      m1=eigen(Sigma)
      v1=sqrt(m1$values)
      vv=diag(v1)
      Pmtx=m1$vectors
      Sh=Pmtx%*%vv%*%t(Pmtx)
   }
   #
   if(k < 1) k=1
   PSI=diag(rep(1,k))
   if(orth){
      WGT=c(PSI%*%Sh)
   }
   else{
      WGT=c(PSI)
   }
   #
   for (il in 1:lag){
      ilk=il*k
      tmp=matrix(0,k,k)
      if((q > 0) && (il <= q))tmp=-Theta[,(ilk-k+1):ilk]
      if(p > 0){
         iend=min(il,p)
         for (j in 1:iend){
            jdx=(il-j)
            kdx=j*k
            tmp=tmp+Phi[,(kdx-k+1):kdx]%*%PSI[,(jdx*k+1):(jdx*k+k)]
         }
         ## end  of p > 0.
      }
      PSI=cbind(PSI,tmp)
      if(orth){
         WGT=cbind(WGT,c(tmp%*%Sh))
      }
      else{
         WGT=cbind(WGT,c(tmp))
      }
      ### end of il-loop
   }
   wk1=WGT
   for (i in 1:k^2){
      wk1[i,] = cumsum(WGT[i,])
   }
   ## plots the psi-weights
   tdx=c(1:(lag+1))-1
   par(mfcol=c(k,k),mai=c(0.3,0.3,0.3,0.3))
   if(orth){
      gmax=max(WGT)
      gmin=min(WGT)
      cx=(gmax-gmin)/10
      gmax=gmax+cx
      gmin=gmin-cx
      for (j in 1:k^2){
         plot(tdx,WGT[j,],type='l',xlab='lag',ylab='IRF',ylim=c(gmin,gmax),cex.axis=0.8)
         points(tdx,WGT[j,],pch='*',cex=0.8)
         title(main='Orth. innovations')
      }
      cat("Press return to continue ","\n")
      readline()
      gmax=max(wk1)
      gmin=min(wk1)
      cx=(gmax-gmin)/10
      gmax=gmax+cx
      gmin=gmin-cx
      for (j in 1:k^2){
         plot(tdx,wk1[j,],type='l',xlab='lag',ylab="Acu-IRF",ylim=c(gmin,gmax),cex.axis=0.8)
         points(tdx,wk1[j,],pch="*",cex=0.8)
         title(main='Orth. innovations')
      }
   }
   else{
      gmax=max(WGT)
      gmin=min(WGT)
      cx=(gmax-gmin)/10
      gmax=gmax+cx
      gmin=gmin-cx
      for(j in 1:k^2){
         plot(tdx,WGT[j,],type='l',xlab='lag',ylab='IRF',ylim=c(gmin,gmax),cex.axis=0.8)
         points(tdx,WGT[j,],pch='*',cex=0.8)
         title(main="Orig. innovations")
      }
      cat("Press return to continue ","\n")
      readline()
      gmax=max(wk1)
      gmin=min(wk1)
      cx=(gmax-gmin)/10
      gmax=gmax+cx
      gmin=gmin-cx
      for(j in 1:k^2){
         plot(tdx,wk1[j,],type='l',xlab='lag',ylab='Acu-IRF',ylim=c(gmin,gmax),cex.axis=0.8)
         points(tdx,wk1[j,],pch='*',cex=0.8)
         title(main="Orig. innovations")
      }
   }
   par(mfcol=c(1,1))
   VARMAirf <- list(psi=PSI,irf=WGT)
 }

"VARMAcov" <- function(Phi=NULL,Theta=NULL,Sigma=NULL,lag=12,trun=120){
   ## trun: trunction point for psi-weights used in the calculation.
   ##
   m1=PSIwgt(Phi=Phi,Theta=Theta,lag=trun,plot=FALSE)
   Psi=m1$psi.weight
   nc=dim(Psi)[2]; k=dim(Psi)[1]
   if(is.null(Sigma)){
      wk=Psi
   }
   else{
      wk=NULL
      for (i in 0:trun){
         ist=i*k
         wk=cbind(wk,Psi[,(ist+1):(ist+k)]%*%Sigma)
      }
      #end of else
   }
   Gam0=wk%*%t(Psi)
   SE=diag(1/sqrt(diag(Gam0)))
   covmtx=Gam0; cormtx=SE%*%Gam0%*%SE
   for (i  in 1:lag){
      ist=i*k
      Gami=wk[,(ist+1):nc]%*%t(Psi[,1:(nc-ist)])
      covmtx=cbind(covmtx,Gami)
      cormtx=cbind(cormtx,SE%*%Gami%*%SE)
   }
   for (i in 0:lag){
      ist=i*k
      cat("Auto-Covariance matrix of lag: ",i,"\n")
      print(round(covmtx[,(ist+1):(ist+k)],5))
   }
   for (i in 0:lag){
      ist=i*k
      cat("cross correlation matrix of lag: ",i,"\n")
      print(round(cormtx[,(ist+1):(ist+k)],4))
   }
   
   VARMAcov <- list(autocov=covmtx,ccm=cormtx)
}

"Eccm" <- function(zt,maxp=5,maxq=6,include.mean=FALSE,rev=TRUE){
   ### Compute extended cross-correlation matrices using iterated regression
   ### fitting instead of the recursive method.
   #### rev: a switch to compute Q(m) statistics from q to maxq.
   if(!is.matrix(zt))zt=as.matrix(zt)
   x=zt
   nT=dim(x)[1]
   k=dim(x)[2]
   if(include.mean){
      av=apply(x,2,mean)
      for (i in 1:k){
         x[,i]=x[,i]-av[i]
      }
   }
   ksq=k*k
   if(rev){
      maxq1=maxq+1
      m1=revmq(x,maxq1,output=F)
      vEccm=m1$ccm[,1:maxq1]
      pEccm=m1$pvalue
      ARcoef=NULL
      for (p in 1:maxp){
         Phi=NULL
         m1=VAR(x,p,include.mean=F,output=F)
         Phi=rbind(Phi,m1$Phi)
         resi=m1$residuals
         m2=revmq(resi,maxq1,output=F)
         pv1=m2$pvalue[1]
         Eccmit=m2$ccm[,1]
         y=x[(p+1):nT,]
         xreg=NULL
         for (j in 1:p){
            xreg=cbind(xreg,x[(p+1-j):(nT-j),])
            }
         kx=dim(xreg)[2]
         for (it in 1:maxq){
            y=y[-1,]
            yT=dim(y)[1]
            if(it == 1){
               xreg=xreg[-1,]
            }
            else{
               nx=dim(xreg)[2]
               xT=dim(xreg)[1]
               xreg=cbind(xreg[-1,1:kx],xreg[-xT,(kx+1):nx])
            }
            TT=dim(resi)[1]
            xreg=cbind(xreg,resi[-TT,])
            xpx=crossprod(xreg,xreg)/nT
            xpy=crossprod(xreg,y)/nT
            xpxinv=solve(xpx)
            beta=xpxinv%*%xpy
            ##cat("beta","\n")
            ##print(beta)
            wt=y-xreg[,1:kx]%*%beta[1:kx,]
            resi=y-xreg%*%beta
            Phi=rbind(Phi,t(beta[1:kx,]))
            m3=revmq(wt,maxq1,output=F)
            Eccmit=cbind(Eccmit,m3$ccm[,(it+1)])
            pv1=c(pv1,m3$pvalue[it+1])
            ##end of it-loop
         }
         ARcoef=cbind(ARcoef,Phi)
         vEccm=rbind(vEccm,Eccmit)
         pEccm=rbind(pEccm,pv1)
       }
    }
   else {
       m1=ccm(x,(maxq+1),output=F)
      vEccm=m1$ccm[,2:(maxq+2)]
      pEccm=m1$pvalue
       ARcoef=NULL
      for (p in 1:maxp){
         Phi=NULL
         m1=VAR(x,p,include.mean=F,output=F)
         Phi=rbind(Phi,m1$Phi)
         resi=m1$residuals
         m2=ccm(resi,1,output=F)
         pv1=m2$pvalue
         Eccmit=matrix(m2$ccm[,2],ksq,1)
         y=x[(p+1):nT,]
         xreg=NULL
         for (j in 1:p){
            xreg=cbind(xreg,x[(p+1-j):(nT-j),])
         }
          kx=dim(xreg)[2]
         for (it in 1:maxq){
            y=y[-1,]
            yT=dim(y)[1]
            if(it == 1){
               xreg=xreg[-1,]
            }
            else{
               nx=dim(xreg)[2]
               xT=dim(xreg)[1]
               xreg=cbind(xreg[-1,1:kx],xreg[-xT,(kx+1):nx])
            }
            TT=dim(resi)[1]
            xreg=cbind(xreg,resi[-TT,])
            xpx=crossprod(xreg,xreg)/nT
            xpy=crossprod(xreg,y)/nT
            xpxinv=solve(xpx)
            beta=xpxinv%*%xpy
             wt=y-xreg[,1:kx]%*%beta[1:kx,]
            resi=y-xreg%*%beta
            Phi=rbind(Phi,t(beta[1:kx,]))
            m3=ccm(wt,it+1,output=F)
            Eccmit=cbind(Eccmit,m3$ccm[,(it+2)])
            pv1=c(pv1,m3$pvalue[it+1])
            ##end of it-loop
         }
          ARcoef=cbind(ARcoef,Phi)
         vEccm=rbind(vEccm,Eccmit)
         pEccm=rbind(pEccm,pv1)
       }
    }
   cat("p-values table of Extended Cross-correlation Matrices:","\n")
   cat("Column: MA order","\n")
   cat("Row   : AR order","\n")
   colnames(pEccm) <- c(c(0:maxq))
   rownames(pEccm) <- c(c(0:maxp))
   tmp=round(pEccm,4)
   printCoefmat(tmp)
    Eccm <- list(pEccm=pEccm,vEccm=vEccm,ARcoef=ARcoef)
  }


"revmq" <- function(x,lag=12,output=FALSE){
   # Compute multivariate Ljung-Box test statistics
   ## Show the results based on the following reversed test:
   ## H_0: rho_{i} = rho_{i+1} = ... = rho_{maxq} = 0.
   ## for i = 1, 2, ...., maxq.
   if(!is.matrix(x))x=as.matrix(x)
   nr=dim(x)[1]
   nc=dim(x)[2]
   nr1=nr-1
   nrsq=nr*nr
   ksq=nc*nc
   x1=scale(x,center=TRUE,scale=FALSE)
   g0=crossprod(x1,x1)/nr1
   S1=sqrt(diag(g0))
   D=diag(1/S1)
   ginv=solve(g0)
   Qm=NULL
   ccm = NULL
   qm=0.0
   for (i in 1:lag){
      x1a=x1[(i+1):nr,]
      x2a=x1[1:(nr-i),]
      g = crossprod(x1a,x2a)/nr1
      rho=D%*%g%*%D
      ccm=cbind(ccm,matrix(c(rho),ksq,1))
      h=t(g)%*%ginv%*%g%*%ginv
      qm=qm+nrsq*sum(diag(h))/(nr-i)
      Qm=c(Qm,qm)
   }
   df=ksq*lag
   rqm=Qm[lag]
   pvs=1-pchisq(Qm[lag],df)
   if(lag > 1){
      for (i in 1:(lag-1)){
         tst=Qm[lag]-Qm[i]
         df=df-ksq
         rqm=c(rqm,tst)
         pvs=c(pvs,1-pchisq(tst,df))
      }
   }
   if(output){
      cat("Qm:","\n")
      print(Qm)
      cat("reversed-qm","\n")
      print(rqm)
      cat("p-values of rqm: ","\n")
      print(pvs)
   }
   revmq <- list(ccm=ccm,rqm=rqm,pvalue=pvs)
 }

"Kronid" <- function(x,plag=5,crit=0.05){
   # Identifies the Kronecker indexes for a vector time series
   # plag is the number of lags used to represent the PAST vector
   if(!is.matrix(x))x=as.matrix(x)
   nT=dim(x)[1]
   k=dim(x)[2]
   y=as.matrix(x)
   if(plag < 1){
      plag=floor(log(nT))+1
   }
   # Construct the PAST-vector
   iend=nT-plag
   past=y[1:iend,]
   if (plag > 1){
      for (i in 2:plag){
         past=cbind(y[i:(iend+i-1),],past)
      }
   }
   # initialize the Kronecker indexes and control variable.
   kdx=rep(0,k)
   found=rep(0,k)
   h=0
   ist=plag+1
   futu1=as.matrix(y[ist:nT,])
   cat("h = ",h,"\n")
   #print(h)
   for (i in 1:k){
      cat("Component = ",i,"\n")
      s1=c(i)
      if(i > 1){
         fnd=found[1:(i-1)]
         jdx=c(1:(i-1))[fnd==0]
         s1=c(jdx,i)
      }
      futu=as.matrix(futu1[,s1])
      m1=cancor(past,futu)
      df=dim(futu)[2]
      dp=dim(past)[2]
      deg=dp-df+1
      seig=m1$cor[df]^2
      cat("square of the smallest can. corr. = ",seig,"\n")
      tst=-(nT-1-0.5*(dp+df-1))*log(1-seig)
      pv=1-pchisq(tst,deg)
      stat=c(tst,deg,pv)
      cat("    test,   df, &  p-value:","\n")
      print(round(stat,3))
      if(i>1){
         cstar=cbind(cstar,stat)
      }
      else{
         cstar=stat
      }
      if(pv > crit){
         found[i]=1
         kdx[i]=h
         cat("A Kronecker index found","\n")
      }
   }
   cat("=============","\n")
   while(sum(found) < k){
      idim=dim(past)[1]
      h=h+1
      cat("h = ",h,"\n")
      past=past[1:(idim-1),]
      futu=futu[1:(idim-1),]
      futu1=y[(ist+h):nT,]
      for (ii in 1:k){
         if(found[ii]==0){
            cat("Component = ",ii,"\n")
            #print(ii)
            futu=cbind(futu,futu1[,ii])
            m1=cancor(past,futu)
            df=dim(futu)[2]
            dp=dim(past)[2]
            deg=dp-df+1
            seig=m1$cor[df]^2
            cat("Square of the smallest can. corr. = ",seig,"\n")
            y1=futu%*%(m1$ycoef[,df])
            x1=past%*%(m1$xcoef[,df])
            m2=acf(y1,lag.max=h,plot=F)
            acfy=m2$acf[2:(h+1)]
            m3=acf(x1,lag.max=h,plot=F)
            acfx=m3$acf[2:(h+1)]
            dsq=1+2*sum(acfx*acfy)
            seig=seig/dsq
            tst=-(nT-1-0.5*(dp+df-1))*log(1-seig)
            pv=1-pchisq(tst,deg)
            stat=c(tst,deg,pv,dsq)
            cat("    test,     df, p-value & d-hat:","\n")
            print(round(stat,3))
            stat=stat[1:3]
            cstar=cbind(cstar,stat)
            if(pv > crit){
               found[ii]=1
               kdx[ii]=h
               futu=futu[,1:(df-1)]
               cat("A Kronecker found","\n")
            }
         }
      }
      cat("============","\n")
   }
   cat("   ","\n")
   cat("Kronecker indexes identified:","\n")
   print(kdx)
   Kronid<-list(index=kdx,tests=cstar)
 }

"Kronspec" <- function(kdx,output=TRUE){
   # Specify a VARMA model for a given set of Kronecker indices
   #
   # Output: 1 = No estimation (the coefficient of Z_{it}
   #         2 = estimation
   #         0 = fixed to zero.
   k=length(kdx)
   if(output){
      cat("Kronecker indices: ",kdx,"\n")
      cat("Dimension: ",k,"\n")
   }
   KK=sort(kdx,index.return=T)
   idx=KK$ix
   p=KK$x[k]
   q=p
   mx=(p+1)*k
   Theta=matrix(2,k,k*(q+1))
   for (i in 1:k){
      Theta[i,i]=1
      if(kdx[i] < q){
         jj=(kdx[i]+1)*k
         Theta[i,(jj+1):mx]=0
      }
   }
   if(k > 1){
      for (i in 1:(k-1)){
         Theta[i,(i+1):k]=0
      }
   }
   # Indicator matrix of AR polynomial: 1 = estimation, 0 denotes zero.
   Phi=Theta
   # specify the Phi(0) lower triangular part
   if (k > 1){
      for (i in 2:k){
         for (j in 1:(i-1)){
            if(kdx[j] <= kdx[i]) Phi[i,j]=0
            # j-loop
         }
         # i-loop
      }
      # for the case of k > 1
   }
   Theta[1:k,1:k]=Phi[1:k,1:k]
   ## specify redundant parameters
   for (i in 1:k){
      for (j in 1:k){
         if(kdx[i] > kdx[j]){
            for (ii in 1:(kdx[i]-kdx[j]))
            Phi[i,ii*k+j]=0
         }
      }
   }
   if(output){
      cat("Notation: ","\n")
      cat(" 0: fixed to 0","\n")
      cat(" 1: fixed to 1","\n")
      cat(" 2: estimation","\n")
      cat("AR coefficient matrices: ","\n")
      print(Phi)
      cat("MA coefficient matrices: ","\n")
      print(Theta)
   }
   Kronspec <- list(PhiID=Phi,ThetaID=Theta)
 }

"Kronfit" <- function(da,kidx,include.mean=T,fixed=NULL,Kpar=NULL,seKpar=NULL,prelim=F,details=F,thres=1.0){
   # Estimation of a vector ARMA model using conditional MLE (Gaussian dist)
   #  The model is specified via Kronecker indices.
   # When prelim=TRUE, fixed is assigned based on the results of AR approximation.
   if(!is.matrix(da))da=as.matrix(da)
   nT=dim(da)[1]; k=dim(da)[2]
    k1=length(kidx)
      if(k1 <= 0){
         k1=k; kidx=rep(1,k)
       }
      maxk=max(kidx)
      m0=Kronspec(kidx,output=F)
      ARid=m0$PhiID; MAid=m0$ThetaID
      #
      print(ARid)
      print(MAid)

 iniKro <- function(da,at,ARid,MAid,include.mean){
   #### z(t) = xi0 z(t) + SUM[xii *z(t-i)] + SUM[Omegai*a(t-i)] + a(t).
   if(!is.matrix(da))da=as.matrix(da)
   if(!is.matrix(at))at=as.matrix(at)
   nT=dim(da)[1];   k=dim(da)[2]
   ## obtain the maximum index value.
   p=floor(dim(ARid)[2]/k)-1
   if(p <= 0)p=1
   ist = p + 1
   ## est: stores the estimates (equation 1, equation 2, etc.)
   est=NULL
   estse=NULL
   for (i in 1:k){
      X=NULL
      Y=da[ist:nT,i]
      if(include.mean)X=rep(1,(nT-p))
      ### This is for i > 1 Only.
      if(i > 1){
         for (j in 1:(i-1)){
            if(ARid[i,j] > 1){
               tmp=at[ist:nT,j]-da[ist:nT,j]
               X=cbind(X,tmp)
            }
         }
       }
      ### setup the lagged AR variables
      for(lag in 1:p){
         jst=lag*k
         for (j in 1:k){
            if(ARid[i,jst+j] > 1){
               tmp=da[(ist-lag):(nT-lag),j]
               X=cbind(X,tmp)
             }
          }
        }
       for(lag in 1:p){
         jst=lag*k
         for (j in 1:k){
            if(MAid[i,jst+j] > 1){
               tmp=at[(ist-lag):(nT-lag),j]
               X=cbind(X,tmp)
             }
          }
       }
      XPX=crossprod(X,X)/nT
      XPXinv=solve(XPX)
      XPY=crossprod(X,Y)/nT
      beta=XPXinv%*%XPY
      l1=dim(XPX)[1]
      resi=Y-X%*%matrix(beta,l1,1)
      evar=crossprod(resi,resi)/(nT-p)
      est=c(est,beta)
      estse=c(estse,sqrt(diag(XPXinv)*evar/nT))
     }
   iniKro <- list(par=est,se=estse)
  }
 
  if(length(Kpar) < 1){
      m1=VARorder(da,maxk+9,output=FALSE)
      porder=m1$aicor
      if(porder < 1)porder=1
      m2=VAR(da,porder,output=FALSE)
      y=da[(porder+1):nT,]
      x=m2$residuals
      m3=iniKro(y,x,ARid,MAid,include.mean)
      ### Kpar is the vector of ALL estimable parameters.
      Kpar <- m3$par; seKpar=m3$se
      ### Kpar is a vector; which stores parameters equation-by-equation.
      nr=length(Kpar) 
      ### Preliminary simplification
      if(prelim){
         fixed = rep(0,nr)
         for (j in 1:nr){
            tt=Kpar[j]/seKpar[j]
            if(abs(tt) >= thres){
               fixed[j]=1
            }
            else{
               Kpar[j]=0
            }
          }
        }
      }
     else{
      nr=length(Kpar)
     }
   if(length(fixed) < 1){fixed=rep(1,nr)}
   # Identify parameters to be estimated.
   nr=length(Kpar)
   JJdx=c(1:nr)[fixed==1]
   par=Kpar[JJdx]
   separ= seKpar[JJdx]
   #########
   cat("Number of parameters: ",length(par),"\n")
   cat("initial estimates: ",round(par,4),"\n")
   ### Set up lower and upper bounds
   lowerBounds=par; upperBounds=par
   for (j in 1:length(par)){
      lowerBounds[j] = par[j]-2*separ[j]
      upperBounds[j] = par[j]+2*separ[j]
   }
   cat("Upper-bound: ",round(upperBounds,4),"\n")
   cat("Lower-bound: ",round(lowerBounds,4),"\n")
### likelihood function 
 LLKron <- function(par,zt=da,JJdx=JJdx,kidx=kidx,ARid=ARid,MAid=MAid,Kpar=Kpar,include.mean=include.mean){
   k=dim(zt)[2]
   nT=dim(zt)[1]
   maxk = max(kidx)
   Kpar[JJdx]=par
   ###  Assign parameters to their proper locations in the program.
   Cnt=rep(0,k)
   Ph0=diag(rep(1,k))
   kp1= dim(ARid)[2]; kp=kp1-k
   PH=matrix(0,k,kp)
   TH=matrix(0,k,kp)
   icnt=0
   for (i in 1:k){
      idx=c(1:kp1)[ARid[i,] > 1]; jdx=c(1:kp1)[MAid[i,] > 1]
      # kdx denotes the number of non-zero elements in lag-0.
      kdx=c(1:k)[ARid[i,1:k] > 1]
      if(length(kdx) > 0){
         idx=idx[-kdx]; jdx=jdx[-kdx]
      }
      iend=length(idx); jend=length(jdx); kend=length(kdx)
      #### icnt: parameter count
      if(include.mean){
         icnt=icnt+1
         Cnt[i]=Kpar[icnt]
      }
      if(kend > 0){
         Ph0[i,kdx]=Kpar[(icnt+1):(icnt+kend)]
         icnt=icnt+kend
      }
      if(iend > 0){
         PH[i,idx-k]=Kpar[(icnt+1):(icnt+iend)]
         icnt=icnt+iend
      }
      if(jend > 0){
         TH[i,jdx-k]=Kpar[(icnt+1):(icnt+jend)]
         icnt=icnt+jend
      }
   }
   Ph0i=solve(Ph0)
   ARc=Ph0i%*%PH
   MAc=Ph0i%*%TH
   Cntc=Ph0i%*%as.matrix(Cnt,k,1)
   ist=maxk+1
   at=matrix((zt[1,]-Cntc),1,k)
   if(maxk > 1){
      for (t in 2:maxk){
         tmp=matrix((zt[t,]-Cntc),1,k)
         for (j in 1:maxk){
            if((t-j) > 0){
               jdx=(j-1)*k
               tmp1=matrix(zt[(t-j),],1,k)%*%t(as.matrix(ARc[,(jdx+1):(jdx+k)]))
               tmp=tmp-tmp1
            }
         }
         for (j in 1:maxk){
            jdx=(j-1)*k
            if((t-j)>0){
               tmp2=matrix(at[(t-j),],1,k)%*%t(as.matrix(MAc[,(jdx+1):(jdx+k)]))
               tmp=tmp-tmp2
            }
           }
         at=rbind(at,tmp)
        }
     }
   ### for t from ist on
   ist=maxk+1
   Pcnt = NULL; beta=NULL
   if(include.mean)beta=matrix(Cntc,1,k)
   beta=rbind(beta,t(ARc),t(MAc))
   idim=k*maxk*2
   if(include.mean){
      Pcnt=c(1)
      idim=idim+1
   }
   #
   for (t in (maxk+1):nT){
      Past=NULL
      for (j in 1:maxk){
         Past=c(Past,zt[(t-j),])
      }
      for (j in 1:maxk){
         Past=c(Past,at[(t-j),])
      }
      tmp = matrix(c(Pcnt,Past),1,idim)%*%beta
      tmp3=zt[t,]-tmp
      at=rbind(at,tmp3)
   }
   at=at[(ist:nT),]
   sig=crossprod(at,at)/(nT-maxk)
   ll=dmvnorm(at,rep(0,k),sig)
   LLKron=-sum(log(ll))
   LLKron
  }
   # Step 5: Estimate Parameters and Compute Numerically Hessian:
   if(details){
      fit = nlminb(start = par, objective = LLKron,zt=da,include.mean=include.mean,JJdx=JJdx,kidx=kidx, 
      Kpar=Kpar,ARid=ARid,MAid=MAid,lower = lowerBounds, upper = upperBounds, control = list(trace=3))
   }
   else {
      fit = nlminb(start = par, objective = LLKron, zt=da,include.mean=include.mean,JJdx=JJdx,kidx=kidx, 
       Kpar=Kpar,ARid=ARid,MAid=MAid,lower = lowerBounds, upper = upperBounds)
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
         Hessian[i, j] = (LLKron(x1,zt=da,include.mean=include.mean,JJdx=JJdx,kidx=kidx,Kpar=Kpar,ARid=ARid,MAid=MAid)
                         -LLKron(x2,zt=da,include.mean=include.mean,JJdx=JJdx,kidx=kidx,Kpar=Kpar,ARid=ARid,MAid=MAid)
                         -LLKron(x3,zt=da,include.mean=include.mean,JJdx=JJdx,kidx=kidx,Kpar=Kpar,ARid=ARid,MAid=MAid)
                         +LLKron(x4,zt=da,include.mean=include.mean,JJdx=JJdx,kidx=kidx,Kpar=Kpar,ARid=ARid,MAid=MAid))/
         (4*epsilon[i]*epsilon[j])
      }
   }
   # Step 6: Create and Print Summary Report:
   d1=det(Hessian)
   if(d1 < 1.0e-10){
      se.coef=rep(1,npar)
   }
   else{
      se.coef = sqrt(diag(solve(Hessian)))
   }
   tval = fit$par/se.coef
   matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
   dimnames(matcoef) = list(names(tval), c(" Estimate",
   " Std. Error", " t value", "Pr(>|t|)"))
   cat("\nCoefficient(s):\n")
   printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
   Kpar[JJdx]=fit$par
   seKpar[JJdx]=se.coef
   # Restore estimates to the format of unconstrained case for printing.
   Cnt=rep(0,k); seCnt=rep(0,k)
   Ph0=diag(rep(1,k)); sePh0=diag(rep(1,k))
   kp1= dim(ARid)[2]; kp=kp1-k
   PH=matrix(0,k,kp); sePH=matrix(0,k,kp)
   TH=matrix(0,k,kp); seTH=matrix(0,k,kp)
   icnt=0
   for (i in 1:k){
      idx=c(1:kp1)[ARid[i,] > 1]; jdx=c(1:kp1)[MAid[i,] > 1]
      # kdx denotes the number of non-zero elements in lag-0.
      kdx=c(1:k)[ARid[i,1:k] > 1]
      if(length(kdx) > 0){
         idx=idx[-kdx]; jdx=jdx[-kdx]
      }
      iend=length(idx); jend=length(jdx); kend=length(kdx)
      if(include.mean){
         icnt=icnt+1
         Cnt[i]=Kpar[icnt]
         seCnt[i]=seKpar[icnt]
      }
      if(kend > 0){
         Ph0[i,kdx]=Kpar[(icnt+1):(icnt+kend)]
         sePh0[i,kdx]=seKpar[(icnt+1):(icnt+kend)]
         icnt=icnt+kend
      }
      if(iend > 0){
         ##cat("idx-k: ",idx-k,"\n")
         PH[i,idx-k]=Kpar[(icnt+1):(icnt+iend)]
         sePH[i,idx-k]=seKpar[(icnt+1):(icnt+iend)]
         icnt=icnt+iend
      }
      if(jend > 0){
         TH[i,jdx-k]=Kpar[(icnt+1):(icnt+jend)]
         seTH[i,jdx-k]=seKpar[(icnt+1):(icnt+jend)]
         icnt=icnt+jend
      }
   }
   cat("---","\n")
   cat("Estimates in matrix form:","\n")
   if(include.mean){
      cat("Constant term: ","\n")
      cat("Estimates: ",round(Cnt,3),"\n")
   }
   cat("AR and MA lag-0 coefficient matrix","\n")
   print(round(Ph0,3))
   cat("AR coefficient matrix","\n")
   jcnt=0
   for (i in 1:maxk){
      cat("AR(",i,")-matrix","\n")
      ph=PH[,(jcnt+1):(jcnt+k)]
      print(round(ph,3))
      jcnt=jcnt+k
   }
   cat("MA coefficient matrix","\n")
   icnt=0
   for (i in 1:maxk){
      cat("MA(",i,")-matrix","\n")
      theta=-TH[,(icnt+1):(icnt+k)]
      print(round(theta,3))
      icnt=icnt+k
   }
   ##### Compute the residuals
   Ph0i=solve(Ph0)
   ARc=Ph0i%*%PH
   MAc=Ph0i%*%TH
   Cntc=Ph0i%*%as.matrix(Cnt,k,1)
   zt=da 
   ist=maxk+1
   #### consider the case t from 1 to maxk+1
   at=matrix((zt[1,]-Cntc),1,k)
   if(maxk > 1){
      for (t in 2:maxk){
         tmp=matrix((zt[t,]-Cntc),1,k)
         for (j in 1:maxk){
            if((t-j) > 0){
               jdx=(j-1)*k
               tmp1=matrix(zt[(t-j),],1,k)%*%t(as.matrix(ARc[,(jdx+1):(jdx+k)]))
               tmp=tmp-tmp1
            }
          }
          for (j in 1:maxk){
            jdx=(j-1)*k
            if((t-j)>0){
               tmp2=matrix(at[(t-j),],1,k)%*%t(as.matrix(MAc[,(jdx+1):(jdx+k)]))
               tmp=tmp-tmp2
            }
          }
         at=rbind(at,tmp)
       }
    }
    ### for t from ist on
   ist=maxk+1
   Pcnt=NULL
   beta=NULL
   if(include.mean){
      beta=matrix(Cntc,1,k)
      Pcnt=c(1)
   }
   beta=rbind(beta,t(ARc),t(MAc))
   idim=k*maxk*2
   if(include.mean){
      Pcnt=c(1)
      idim=idim+1
   }
   #
   for (t in (maxk+1):nT){
      Past=NULL
      for (j in 1:maxk){
         Past=c(Past,zt[(t-j),])
      }
      for (j in 1:maxk){
         Past=c(Past,at[(t-j),])
      }
      tmp = matrix(c(Pcnt,Past),1,idim)%*%beta
      tmp3=zt[t,]-tmp
      at=rbind(at,tmp3)
   }
    at=at[(ist:nT),]
   sig=crossprod(at,at)/(nT-maxk)
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
    
   Kronfit <- list(data=da,Kindex=kidx,ARid=ARid,MAid=MAid,cnst=include.mean,coef=Kpar,secoef=seKpar,residuals=at,Sigma=sig,aic=aic,bic=bic, Ph0=Ph0,Phi=PH,Theta=-TH)
}

"refKronfit" <- function(model,thres=1.0){
   zt=model$data
   inc.mean=model$cnst
   kidx=model$Kindex
   Kpar= model$coef
   seKpar= model$secoef
   maxk=max(kidx)
   nr=length(Kpar)
   fix=rep(0,nr)
   for (j in 1:nr){
      tt = 0
      iav=is.na(seKpar[j])
      if(iav)seKpar[j]=0.01
      tt=Kpar[j]/seKpar[j]
      if(abs(tt) > thres){
         fix[j]=1
      }
      else{
         Kpar[j]=0
      }
     }
   m1=Kronfit(zt,kidx,include.mean=inc.mean,fixed=fix,Kpar=Kpar,seKpar=seKpar)
   ARid=m1$ARid; MAid=m1$MAid
   Kpar=m1$coef; seKpar=m1$secoef
   sig=m1$Sigma; aic=m1$aic; bic=m1$bic
   Ph0=m1$Ph0
   PH=m1$Phi
   TH=-m1$Theta
   at=m1$residuals
   
   refKronfit <- list(data=zt,Kindex=kidx,ARid=ARid,MAid=MAid,cnst=inc.mean,coef=Kpar,secoef=seKpar,residuals=at,Sigma=sig,aic=aic,bic=bic, Ph0=Ph0,Phi=PH,Theta=-TH)
}

"sVARMA" <- function(da,order=c(0,0,0),sorder=c(0,0,0),s=12,include.mean=T,fixed=NULL,details=F,switch=F){
   # Estimation of a multiplicative vector ARMA model using conditional MLE (Gaussian dist)
   if(!is.matrix(da))da=as.matrix(da)
   p=order[1];d=order[2];q=order[3];P=sorder[1];D=sorder[2];Q=sorder[3]
   nT=dim(da)[1]; k=dim(da)[2]
   # basic setup.
   if(p < 0)p=0; if(q < 0)q=0; if(P < 0) P = 0; if(Q < 0) Q = 0; if(s < 0) s=-s
   if(d > 1){
      cat("Regular difference is adjusted to d=1","\n")
      d=1
   }
   if(D > 1){
      cat("Seasonal difference is adjusted to D=1","\n")
      D=1
   }
   kp=k*p
   kq=k*q
   kP=k*P
   kQ=k*Q
   # Take care of the difference
   if(d==1){
      X=NULL
      MEAN=rep(0,k)
      for (j in 1:k){
         X=cbind(X,diff(da[,j]))
         t1=t.test(X[,j])
         if(t1$p.value < 0.05)MEAN[j]=1
      }
      if(sum(MEAN) < 1)include.mean=FALSE
   }
   else{
      X=da
   }
   if(D==1){
      DX=NULL
      Smean=rep(0,k)
      for (j in 1:k){
         DX=cbind(DX,diff(X[,j],s))
         t1=t.test(DX[,j])
         if(t1$p.value < 0.05)Smean[j]=1
      }
      if(sum(Smean) < 1)include.mean=FALSE 
   }
   else{
      DX=X
   }
   nT=dim(DX)[1]
   arlags=NULL
   if(p > 0){
      arlags=c(1:p)
      if(P > 0)arlags=c(arlags,c(1:P)*s,c(1:P)*s+c(1:p))
   }
   else{
      if(P > 0)arlags=c(1:P)*s
   }
   malags=NULL
   if(q > 0){
      malags=c(1:q)
      if(Q > 0)malags=c(malags,c(1:Q)*s,c(1:Q)*s+c(1:q))
   }
   else{
      if(Q > 0)malags=c(1:Q)*s
   }
   # number of AR and MA lags of the model
   nar=length(arlags)
   nma=length(malags)
   idim=k*(nar+nma)
   if(include.mean)idim=idim+1
   if(length(fixed)==0){fixed=matrix(1,idim,k)}
   Order <- c(order,sorder)
   ARlags <- arlags; MAlags <- malags
   ####
   phi=NULL; sphi=NULL; sephi=NULL; sesphi=NULL
   if(p > 0)phi=matrix(0,k,k*p); sephi=phi
   if(P > 0)sphi=matrix(0,k,k*P);sesphi=sphi
   theta=NULL; stheta=NULL;setheta=NULL; sestheta=NULL
   if(q > 0)theta=matrix(0,k,k*q);setheta=theta
   if(Q > 0)stheta=matrix(0,k,k*Q);sestheta=stheta
   ## Obtain initial estimates of the component parameters using univariate models.
   ### For cross-series initial estimates, we use linear models with univariate at-series
   resi=NULL
   for (j in 1:k){
      m1=arima(DX[,j],order=c(p,0,q),seasonal=list(order=c(P,0,Q),period=s))
      resi=cbind(resi,m1$residuals)
      seest=sqrt(diag(m1$var.coef))
      icnt=0
      if(p > 0){
         for (i in 1:p){
            icnt=icnt+1
            ii=(i-1)*k
            phi[j,(ii+j)]=m1$coef[icnt]
            sephi[j,(ii+j)]=seest[icnt]
         }
      }
      if(q > 0){
         for (i in 1:q){
            ii=(i-1)*k
            icnt=icnt+1
            theta[j,(ii+j)]=-m1$coef[icnt]
            setheta[j,(ii+j)]=seest[icnt]
         }
      }
      if(P > 0){
         for (i in 1:P){
            icnt=icnt+1
            ii=(i-1)*k
            sphi[j,(ii+j)]=m1$coef[icnt]
            sesphi[j,(ii+j)]=seest[icnt]
         }
      }
      if(Q > 0){
         for (i in 1:Q){
            ii=(i-1)*k
            icnt=icnt+1
            stheta[j,(ii+j)]=-m1$coef[icnt]
            sestheta[j,(ii+j)]=seest[icnt]
         }
      }
    }

siniEST <- function(y,x,arlags,malags,include.mean){
   if(!is.matrix(y))y=as.matrix(y)
   if(!is.matrix(x))x=as.matrix(x)
   nT=dim(y)[1]
   k=dim(y)[2]
   nar=length(arlags)
   nma=length(malags)
   p=0; if(nar > 0)p=arlags[nar]
   q=0; if(nma > 0)q=malags[nma]
   pq=max(p,q)
   ist=1+pq
   ne=nT-pq
   if(include.mean){
      xmtx=matrix(1,ne,1)
   }
   else {
      xmtx=NULL
   }
   ymtx=as.matrix(y[ist:nT,])
   if(nar > 0){
      for (j in 1:nar){
         jj=arlags[j]
         xmtx=cbind(xmtx,y[(ist-jj):(nT-jj),])
      }
   }
   if(nma > 0){
      for (j in 1:nma){
         jj=malags[j]
         xmtx=cbind(xmtx,x[(ist-jj):(nT-jj),])
      }
   }
   xmtx=as.matrix(xmtx)
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
   siniEST <- list(estimates=beta,se=sebeta)
 }
   #### Obtain estimates of cross-series parameters, using Least-Squares approximation.
   m2=siniEST(DX,resi,arlags,malags,include.mean)
   #### Fill in the coefficient matrices
   beta=t(m2$estimates)
   sebeta=t(m2$se)
   ##
   icnst=0
   if(include.mean)icnst=1
   if(nar > 0){
      if(p > 0){
         for (i in 1:p){
            idx=(i-1)*k
            for (ii in 1:k){
               jdx=c(1:k)[-ii]
               phi[jdx,(idx+ii)]=beta[jdx,(icnst+idx+ii)]
               sephi[jdx,(idx+ii)]=sebeta[jdx,(icnst+idx+ii)]
            }
         }
       }
      if(P > 0){
         for (i in 1:P){
            kdx=(i-1)*k
            idx=k*p+kdx
            for (ii in 1:k){
               jdx=c(1:k)[-ii]
               sphi[jdx,(kdx+ii)]=beta[jdx,(icnst+idx+ii)]
               sesphi[jdx,(kdx+ii)]=sebeta[jdx,(icnst+idx+ii)]
            }
         }
      }
    }
   if(nma > 0){
      if(q > 0){
         for (i in 1:q){
            kdx=(i-1)*k
            idx=nar*k+kdx
            for (ii in 1:k){
               jdx=c(1:k)[-ii]
               theta[jdx,(kdx+ii)]=-beta[jdx,(icnst+idx+ii)]
               setheta[jdx,(kdx+ii)]=sebeta[jdx,(icnst+idx+ii)]
            }
         }
      }
      if(Q > 0){
         for (i in 1:Q){
            kdx=(i-1)*k
            idx=(nar+q)*k+kdx
            for (ii in 1:k){
               jdx=c(1:k)[-ii]
               stheta[jdx,(kdx+ii)]=-beta[jdx,(icnst+idx+ii)]
               sestheta[jdx,(kdx+ii)]=sebeta[jdx,(icnst+idx+ii)]
            }
         }
      }
    }
   # Identify parameters to be estimated.
   par=NULL
   separ=NULL
   ist=0
   ## We took the transpose of beta and sebeta after siniEST program.
   if(include.mean){
      jdx=c(1:k)[fixed[1,]==1]
      if(length(jdx) > 0){
         par=beta[jdx,1]
         separ=sebeta[jdx,1]
      }
      ist=1
    }
   if(nar > 0){
      if(p > 0){
         for (j in 1:k){
            idx=c(1:kp)[fixed[(ist+1):(ist+kp),j]==1]
            if(length(idx) > 0){
               par=c(par,phi[j,idx])
               separ=c(separ,sephi[j,idx])
            }
          }
         ist=ist+kp
      }
     if(P > 0){
         for (j in 1:k){
            idx=c(1:kP)[fixed[(ist+1):(ist+kP),j]==1]
            if(length(idx) > 0){
               par=c(par,sphi[j,idx])
               separ=c(separ,sesphi[j,idx])
            }
         }
         ist=ist+kP
      }
    }
   if(nma > 0){
      if(q > 0){
         for (j in 1:k){
            idx=c(1:kq)[fixed[(ist+1):(ist+kq),j]==1]
            if(length(idx) > 0){
               par=c(par,theta[j,idx])
               separ=c(separ,setheta[j,idx])
            }
         }
         ist=ist+kq
      }
      if(Q > 0){
         for (j in 1:k){
            idx=c(1:kQ)[fixed[(ist+1):(ist+kQ),j]==1]
            if(length(idx) > 0){
               par=c(par,stheta[j,idx])
               separ=c(separ,sestheta[j,idx])
            }
         }
      }
    }
   #### keep the first few residuals to be used in likelihood evaluation to compute "at".
   jst=max(arlags,malags)
   Sresi <- resi[1:jst,]
   cat("Number of parameters: ",length(par),"\n")
   cat("initial estimates: ",par,"\n")
   lowerBounds=par; upperBounds=par
   for (j in 1:length(par)){
      lowerBounds[j] = par[j]-2*separ[j]
      upperBounds[j] = par[j]+2*separ[j]
   }

LLKsvarma <- function(par,zt=DX,Order=Order,ARlags=arlags,MAlags=malags,include.mean=include.mean,fixed=fixed,swi=switch,Sresi=Sresi){
   ## recall the relevant information.
   k <- dim(zt)[2];   nT <- dim(zt)[1]
   p=Order[1];q=Order[3];P=Order[4];Q=Order[6]
   kp=k*p;kP=k*P;kq=k*q;kQ=k*Q
   nar=length(ARlags); nma=length(MAlags)
   istart=max(ARlags,MAlags)+1
   ###  Assign parameters to their proper locations in the program.
   beta=NULL
   ist=0
   icnt=0
   Ph0=rep(0,k)
   if(include.mean){
      idx=c(1:k)[fixed[1,]==1]
      icnt=length(idx)
      if(icnt > 0){
         Ph0[idx]=par[1:icnt]
      }
      ist=1
      beta=rbind(beta,Ph0)
   }
   PH=NULL;sPH=NULL
   if(nar > 0){
      if(p > 0){
         PH = matrix(0,k,kp)
         for (j in 1:k){
            idx=c(1:kp)[fixed[(ist+1):(ist+kp),j]==1]
            jdx=length(idx)
            if(jdx > 0){
               PH[j,idx]=par[(icnt+1):(icnt+jdx)]
               icnt=icnt+jdx
            }
            # end of j-loop
         }
         ist=ist+kp
         #end of if (p > 0)
      }
      #### Seasonal AR part
      if(P > 0){
         sPH=matrix(0,k,kP)
         for (j in 1:k){
            idx=c(1:kP)[fixed[(ist+1):(ist+kP),j]==1]
            jdx=length(idx)
            if(jdx > 0){
               sPH[j,idx]=par[(icnt+1):(icnt+jdx)]
               icnt=icnt+jdx
            }
         }
         ist=ist+kP
       }
    }
   TH=NULL;sTH=NULL
   if(nma > 0){
      if(q > 0){
         TH=matrix(0,k,kq)
         for (j in 1:k){
            idx=c(1:kq)[fixed[(ist+1):(ist+kq),j]==1]
            jdx=length(idx)
            if(jdx > 0){
               TH[j,idx]=par[(icnt+1):(icnt+jdx)]
               icnt=icnt+jdx
            }
          }
         ist=ist+kq
       }
      if(Q > 0){
         sTH=matrix(0,k,kQ)
         for (j in 1:k){
            idx=c(1:kQ)[fixed[(ist+1):(ist+kQ),j]==1]
            jdx=length(idx)
            if(jdx > 0){
               sTH[j,idx]=par[(icnt+1):(icnt+jdx)]
               icnt=icnt+jdx
            }
         }
        }
     }
   # Obtain the product of matrix polynomials if necessary
   if((p > 0)&&(P > 0)){
      if(swi){
         Phi=Mtxprod1(PH,sPH,p,P)
      }
      else{
         Phi=Mtxprod(PH,sPH,p,P)
      }
      beta=rbind(beta,t(Phi))
   }
   if((p > 0)&&(P==0))beta=rbind(beta,t(PH))
   if((p==0)&&(P > 0))beta=rbind(beta,t(sPH))
   #
   if((q > 0)&&(Q > 0)){
      if(swi){
         Theta=Mtxprod1(TH,sTH,q,Q)
      }
      else{
         Theta=Mtxprod(TH,sTH,q,Q)
      }
      beta=rbind(beta,-t(Theta))
   }
   if((q > 0)&&(Q==0))beta=rbind(beta,-t(TH))
   if((q==0)&&(Q > 0))beta=rbind(beta,-t(sTH))
   #
    #### consider the case t from 1 to pqmatx
   at=Sresi
   ### for t from istart to T
   Pcnt = NULL
   idim=k*(nar+nma)
   if(include.mean){
      Pcnt=c(1)
      idim=idim+1
   }
   for (t in istart:nT){
      Past=NULL
      if(nar > 0){
         for (j in 1:nar){
            jj=ARlags[j]
            Past=c(Past,zt[(t-jj),])
         }
      }
      if(nma > 0){
         for (j in 1:nma){
            jj=MAlags[j]
            Past=c(Past,at[(t-jj),])
         }
      }
      tmp = matrix(c(Pcnt,Past),1,idim)%*%beta
      tmp3=zt[t,]-tmp
      at=rbind(at,tmp3)
   }
   at=at[istart:nT,]
   sig=t(at)%*%at/(nT-istart+1)
   ll=dmvnorm(at,rep(0,k),sig)
   LLKsvarma=-sum(log(ll))
####   cat("test: ",LLKsvarma,"\n")
   LLKsvarma
  }
 ## estimation
   if(details){
      fit = nlminb(start = par, objective = LLKsvarma,zt=DX,Order=Order,ARlags=ARlags,MAlags=MAlags,include.mean=include.mean,
       fixed=fixed,swi=switch,Sresi=Sresi,lower = lowerBounds, upper = upperBounds, control = list(trace=3))
   }
   else {
      fit = nlminb(start = par, objective = LLKsvarma,zt=DX,Order=Order,ARlags=ARlags,MAlags=MAlags,include.mean=include.mean, 
          fixed=fixed,swi=switch,Sresi=Sresi,lower = lowerBounds, upper = upperBounds)
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
         Hessian[i, j] = 
          (LLKsvarma(x1,zt=DX,Order=Order,ARlags=ARlags,MAlags=MAlags,include.mean=include.mean,fixed=fixed,swi=switch,Sresi=Sresi)
          -LLKsvarma(x2,zt=DX,Order=Order,ARlags=ARlags,MAlags=MAlags,include.mean=include.mean,fixed=fixed,swi=switch,Sresi=Sresi)
          -LLKsvarma(x3,zt=DX,Order=Order,ARlags=ARlags,MAlags=MAlags,include.mean=include.mean,fixed=fixed,swi=switch,Sresi=Sresi)
          +LLKsvarma(x4,zt=DX,Order=Order,ARlags=ARlags,MAlags=MAlags,include.mean=include.mean,fixed=fixed,swi=switch,Sresi=Sresi))/
         (4*epsilon[i]*epsilon[j])
      }
   }
   # Step 6: Create and Print Summary Report:
   se.coef = sqrt(diag(solve(Hessian)))
   tval = fit$par/se.coef
   matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
   dimnames(matcoef) = list(names(tval), c(" Estimate",
   " Std. Error", " t value", "Pr(>|t|)"))
   cat("\nCoefficient(s):\n")
   printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
   est=fit$par
   ### restore estimates to the format of unconstrained case for printing purpose.
   ist=0
   icnt = 0
   Ph0=rep(0,k)
   sePh0=rep(0,k)
   beta=NULL
   sebeta=NULL
   if(include.mean){
      idx=c(1:k)[fixed[1,]==1]
      icnt=length(idx)
      if(icnt > 0){
         Ph0[idx]=est[1:icnt]
         sePh0[idx]=se.coef[1:icnt]
      }
      ist=1
      beta=rbind(beta,Ph0)
      sebeta=rbind(sebeta,sePh0)
   }
   PH=NULL; sePH=NULL; sPH=NULL; sesPH=NULL
   if(p > 0){
      PH=matrix(0,kp,k)
      sePH=matrix(0,kp,k)
      for (j in 1:k){
         idx=c(1:kp)[fixed[(ist+1):(ist+kp),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            PH[idx,j]=est[(icnt+1):(icnt+jdx)]
            sePH[idx,j]=se.coef[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
         # end of j-loop
      }
      #end of if (p > 0)
      ist=ist+kp
      beta=rbind(beta,PH)
      sebeta=rbind(sebeta,sePH)
   }
   if(P > 0){
      sPH=matrix(0,kP,k)
      sesPH=matrix(0,kP,k)
      for (j in 1:k){
         idx=c(1:kP)[fixed[(ist+1):(ist+kP),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            sPH[idx,j]=est[(icnt+1):(icnt+jdx)]
            sesPH[idx,j]=se.coef[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
      }
       ist=ist+kP
      beta=rbind(beta,sPH)
      sebeta=rbind(sebeta,sesPH)
   }
   TH=NULL;seTH=NULL; sTH=NULL; sesTH=NULL
   if(q > 0){
      TH=matrix(0,kq,k)
      seTH=matrix(0,kq,k)
      for (j in 1:k){
         idx=c(1:kq)[fixed[(ist+1):(ist+kq),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            TH[idx,j]=est[(icnt+1):(icnt+jdx)]
            seTH[idx,j]=se.coef[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
       }
      ist=ist+kq
      beta=rbind(beta,-TH)
      sebeta=rbind(sebeta,seTH)
   }
   if(Q > 0){
      sTH=matrix(0,kQ,k)
      sesTH=matrix(0,kQ,k)
      for (j in 1:k){
         idx=c(1:kQ)[fixed[(ist+1):(ist+kQ),j]==1]
         jdx=length(idx)
         if(jdx > 0){
            sTH[idx,j]=est[(icnt+1):(icnt+jdx)]
            sesTH[idx,j]=se.coef[(icnt+1):(icnt+jdx)]
            icnt=icnt+jdx
         }
      }
      beta=rbind(beta,-sTH)
      sebeta=rbind(sebeta,sesTH)
   }
   cat("---","\n")
   cat("Estimates in matrix form:","\n")
   if(include.mean){
      cat("Constant term: ","\n")
      cat("Estimates: ",Ph0,"\n")
   }
   if(p > 0){
      cat("Regular AR coefficient matrix","\n")
      jcnt=0
      for (i in 1:p){
         cat("AR(",i,")-matrix","\n")
         ph=t(PH[(jcnt+1):(jcnt+k),])
         print(ph,digits=3)
         jcnt=jcnt+k
      }
    }
    if(P > 0){
      cat("Seasonal AR coefficient matrix","\n")
      jcnt=0
      for (i in 1:P){
         cat("AR(",i*s,")-matrix","\n")
         ph=t(sPH[(jcnt+1):(jcnt+k),])
         print(ph,digits=3)
         jcnt=jcnt+k
      }
    }
   if(q > 0){
      cat("Regular MA coefficient matrix","\n")
      icnt=0
      for (i in 1:q){
         cat("MA(",i,")-matrix","\n")
         the=t(TH[(icnt+1):(icnt+k),])
         print(the,digits=3)
         icnt=icnt+k
      }
    }
   if(Q > 0){
      cat("Seasonal MA coefficient matrix","\n")
      icnt=0
      for (i in 1:Q){
         cat("MA(",i*s,")-matrix","\n")
         the=t(sTH[(icnt+1):(icnt+k),])
         print(the,digists=3)
         icnt=icnt+k
      }
    }
   ######### Obtain product coefficient matrices
   if((p > 0)&&(P > 0)){
      if(switch){
         Phi=t(Mtxprod1(t(PH),t(sPH),p,P))
      }
      else{
         Phi=t(Mtxprod(t(PH),t(sPH),p,P))
      }
   }
   if((p > 0)&&(P==0))Phi=PH
   if((p==0)&&(P > 0))Phi=sPH
   if((q > 0)&&(Q > 0)){
      if(switch){
         Theta=t(Mtxprod1(t(TH),t(sTH),q,Q))
      }
      else{
         Theta=t(Mtxprod(t(TH),t(sTH),q,Q))
      }
   }
   #
   if((q > 0)&&(Q==0))Theta=TH
   if((q==0)&&(Q > 0))Theta=sTH
   ##### Compute the residuals
   zt=DX
   pqmax=max(ARlags,MAlags)
   ist=pqmax+1
   #### consider the case t from ist to T
   at=Sresi[1:pqmax,]
   for (t in ist:nT){
      tmp=zt[t,]-Ph0
      if(nar > 0){
         for (j in 1:nar){
            jj=ARlags[j]
            jdx=(j-1)*k
            ph=Phi[(jdx+1):(jdx+k),]
            tmp=tmp-matrix(zt[(t-jj),],1,k)%*%ph
         }
      }
      if(nma > 0){
         for (j in 1:nma){
            jj=MAlags[j]
            jdx=(j-1)*k
            th=Theta[(jdx+1):(jdx+k),]
            tmp=tmp+matrix(at[(t-jj),],1,k)%*%th
         }
      }
      at=rbind(at,tmp)
   }
   at=at[(ist:nT),]
   c1 = rep("resi",k)
   colnames(at) <- c1
   sig=t(at)%*%at/(nT-pqmax)
   cat(" ","\n")
   cat("Residuals cov-matrix:","\n")
   print(sig)
   dd=det(sig)
   d1=log(dd)
   aic=d1+2*npar/nT
   bic=d1+log(nT)*npar/nT
   cat("----","\n")
   cat("aic= ",round(aic,4),"\n")
   cat("bic= ",round(bic,4),"\n")
   if(length(PH) > 0)PH=t(PH)
   if(length(sPH) > 0)sPH=t(sPH)
   if(length(TH) > 0)TH=t(TH)
   if(length(sTH) > 0)sTH=t(sTH)
   
   sVARMA <- list(data=da,order=order,sorder=sorder,period=s,cnst=include.mean,coef=beta,secoef=sebeta,residuals=at,Sigma=sig,aic=aic,bic=bic,regPhi=PH,seaPhi=sPH, regTheta=TH, seaTheta=sTH, Ph0=Ph0,switch=switch)
}

"Mtxprod" <- function(Mtx,sMtx,p,P){
   # obtain the coefficient matrices of product of two matrix polynomials
   if(!is.matrix(Mtx))Mtx=as.matrix(Mtx)
   if(!is.matrix(sMtx))sMtx=as.matrix(sMtx)
   k=dim(Mtx)[1]
   kp=dim(Mtx)[2]
   kP=dim(sMtx)[2]
   #
   pMtx=Mtx
   for (i in 1:P){
      ii=(i-1)*k
      m2=sMtx[,(ii+1):(ii+k)]
      pMtx=cbind(pMtx,m2)
      for (j in 1:p){
         jdx=(j-1)*k
         m1=Mtx[,(jdx+1):(jdx+k)]
         pMtx=cbind(pMtx,-m1%*%m2)
      }
   }
   pMtx
 }

"Mtxprod1" <- function(Mtx,sMtx,p,P){
   # obtain the coefficient matrices of product of two matrix polynomials
   if(!is.matrix(Mtx))Mtx=as.matrix(Mtx)
   if(!is.matrix(sMtx))sMtx=as.matrix(sMtx)
   k=dim(Mtx)[1]
   kp=dim(Mtx)[2]
   kP=dim(sMtx)[2]
   #
   pMtx=Mtx
   for (i in 1:P){
      ii=(i-1)*k
      m2=sMtx[,(ii+1):(ii+k)]
      pMtx=cbind(pMtx,m2)
      for (j in 1:p){
         jdx=(j-1)*k
         m1=Mtx[,(jdx+1):(jdx+k)]
         pMtx=cbind(pMtx,-m2%*%m1)
      }
   }
   pMtx
}

"refsVARMA" <- function(model,thres=0.8){
   # This program refines the fitted models of sVARMA output by removing
   # insigificant parameters with abs(t-ratio) < thres.
   # model: output object from sVARMA
   # thres: threshold value
   #
   x = model$data
   order=model$order
   sorder=model$sorder
   s=model$period
   cnst = model$cnst
   swi = model$switch
   coef=as.matrix(model$coef)
   secoef=as.matrix(model$secoef)
   nr=dim(coef)[1]
   nc=dim(coef)[2]
   for (j in 1:nc){
      idx=is.na(secoef[,j])
      jdx=c(1:nr)[idx==T]
      secoef[jdx,j]=0.01
   }
   fix=matrix(0,nr,nc)
   for (j in 1:nc){
      tt=coef[,j]/secoef[,j]
      idx=c(1:nr)[abs(tt) >= thres]
      fix[idx,j]=1
   }
   ### Try to keep the constant if the t-ratio is greater then 1.
   if(cnst){
      tt=coef[1,]/secoef[1,]
      idx=c(1:nc)[abs(tt) > 1.0]
      if(length(idx) > 0)fix[1,idx]=1
   }
   
   mm=sVARMA(x,order,sorder,s,include.mean=cnst,fixed=fix,switch=swi)
   
   refsVARMA <- list(data=x,coef=mm$coef,secoef=mm$secoef,order=mm$order,sorder=mm$sorder,period=mm$period,cnst=cnst,residuals=mm$residuals,regPhi=mm$regPhi,seaPhi=mm$seaPhi,regTheta=mm$regTheta,seaTheta=mm$seaTheta,Ph0=mm$Pho,Sigma=mm$Sigma,aic=mm$aic,bic=mm$bic,switch=mm$switch)
   
}

"VARX" <- function(zt,p,xt=NULL,m=0,include.mean=T,fixed=NULL,output=T){
   #This command fits the model
   ## z(t) = c0 + sum_{i=1}^p phi_i * z(t-i) + \sum_{j=0}^m xt(t-j) + a(t).
   ##
   zt=as.matrix(zt)
   if(length(xt) < 1){
      m = -1; kx=0}
   else{
      xt=as.matrix(xt); kx=dim(xt)[2]
    }
   if(p < 0)p=0
   ist=max(p,m)+1
   nT=dim(zt)[1]
   k=dim(zt)[2]
   yt=zt[ist:nT,]
   xmtx=NULL
   if(include.mean)xmtx=rep(1,(nT-ist+1))
   #
   if(p > 0){
      for (i in 1:p){
         xmtx=cbind(xmtx,zt[(ist-i):(nT-i),])
      }
   }
   #
   if( m > -1){
      for (j in 0:m){
         xmtx=cbind(xmtx,xt[(ist-j):(nT-j),])
      }
   }
   #
   p1=dim(xmtx)[2]
   nobe=dim(xmtx)[1]
   ##cat("dim of xmtx",c(nobe,p1),"\n")
   #
   if(length(fixed) < 1){
      ## no constriants
      xpx=t(xmtx)%*%xmtx
      xpy=t(xmtx)%*%yt
      xpxi=solve(xpx)
      beta=xpxi%*%xpy
      resi=as.matrix(yt-xmtx%*%beta)
      sig=crossprod(resi,resi)/nobe
      co=kronecker(sig,xpxi)
      se=sqrt(diag(co))
      se.beta=matrix(se,nrow(beta),k)
      npar=nrow(beta)*k
      d1=log(det(sig))
      aic=d1+2*npar/nobe
      bic=d1+(npar*log(nobe))/nobe
   }
   else{
      # with zero-parameter constriants
      beta=matrix(0,p1,k)
      se.beta=matrix(1,p1,k)
      resi=yt
      npar=0
      for (i in 1:k){
         idx=c(1:p1)[fixed[,i] > 0]
         npar=npar+length(idx)
         if(length(idx) > 0){
            xm=as.matrix(xmtx[,idx])
            y1=matrix(yt[,i],nobe,1)
            xpx=t(xm)%*%xm
            xpy=t(xm)%*%y1
            xpxi=solve(xpx)
            beta1=xpxi%*%xpy
            res = y1 - xm%*%beta1
            sig1=sum(res^2)/nobe
            se=sqrt(diag(xpxi)*sig1)
            beta[idx,i]=beta1
            se.beta[idx,i]=se
            resi[,i]=res
         }
         # end of for (i in 1:k)
      }
      #
      sig=crossprod(resi,resi)/nobe
      d1=log(det(sig))
      aic=d1+2*npar/nobe
      bic=d1+log(nobe)*npar/nobe
      # end of else
   }
   ### print
   Ph0=NULL
   icnt=0
   if(include.mean){
      Ph0=beta[1,]; icnt=icnt+1
      cat("constant term: ","\n")
      cat("est: ",round(Ph0,4),"\n")
      cat(" se: ",round(se.beta[1,],4),"\n")
   }
   Phi=NULL
   if(p > 0){
      Phi=t(beta[(icnt+1):(icnt+k*p),])
      sePhi=t(se.beta[(icnt+1):(icnt+k*p),])
      for (j in 1:p){
         cat("AR(",j,") matrix","\n")
         jcnt=(j-1)*k
         print(round(Phi[,(jcnt+1):(jcnt+k)],3))
         cat("standard errors","\n")
         print(round(sePhi[,(jcnt+1):(jcnt+k)],3))
      }
      icnt=icnt+k*p
      ## end of if(p > 0)
   }
   if(m > -1){
      cat("Coefficients of exogenous","\n")
      Beta=t(beta[(icnt+1):(icnt+(m+1)*kx),])
      seBeta=t(se.beta[(icnt+1):(icnt+(m+1)*kx),])
      for (i in 0:m){
         jdx=i*kx
         cat("lag-",i," coefficient matrix","\n")
         print(round(Beta[,(jdx+1):(jdx+kx)],3))
         cat("standard errors","\n")
         print(round(seBeta[,(jdx+1):(jdx+kx)],3))
      }
      ## end of if(m > -1)
   }
   ##
   cat("Residual Covariance Matrix","\n")
   print(round(sig,5))
   cat("===========","\n")
   cat("Information criteria: ","\n")
   cat("AIC: ",aic,"\n")
   cat("BIC: ",bic,"\n")
   
   VARX <- list(data=zt,xt=xt,aror=p,m=m,Ph0=Ph0,Phi=Phi,beta=Beta,residuals=resi,Sigma=sig,
   coef=beta,se.coef=se.beta,include.mean=include.mean)
}

##### Refine VARX model
"refVARX" <- function(m1,thres=1.0){
   zt=m1$data; xt=m1$xt
   p=m1$aror; m=m1$m; include.m=m1$include.mean
   beta=m1$coef; se.beta=m1$se.coef
   fix=matrix(0,nrow(beta),ncol(beta))
   for (i in 1:ncol(beta)){
      tt=beta[,i]/se.beta[,i]
      idx=c(1:nrow(beta))[abs(tt) > thres]
      if(length(idx) > 0){
         fix[idx,i]=1}
   }
   
   mm=VARX(zt,p,xt,m,include.mean=include.m,fixed=fix)
   Ph0=mm$Ph0; Phi=mm$Phi; Beta=mm$beta; resi=mm$residuals
   sig=mm$Sigma; coef=mm$beta; se.coef=mm$se.coef
   
   refVARX <- list(data=zt,aror=p,xt=xt,m=m,Ph0=Ph0,Phi=Phi,beta=Beta,residuals=resi,Sigma=sig,
   coef=beta,se.coef=se.beta,include.mean=include.m)
}

##### Prediction of VARX models.
#####
"VARXpred" <- function(m1,newxt=NULL,hstep=1,orig=0){
   #This program predicts the VARX model.
   ## z(t) = c0 + sum_{i=1}^p phi_i * z(t-i) + \sum_{j=0}^m xt(t-j) + a(t).
   ##
   zt=as.matrix(m1$data); xt=as.matrix(m1$xt); p=m1$aror; m=m1$m
   Ph0=as.matrix(m1$Ph0); Phi=as.matrix(m1$Phi); Sig=as.matrix(m1$Sigma); beta=as.matrix(m1$beta)
   include.m=m1$include.mean
   nT=dim(zt)[1]; k=dim(zt)[2]; dx=dim(xt)[2]
   se=NULL
   if(length(Ph0) < 1)Ph0=matrix(rep(0,k),k,1)
   if(hstep < 1)hstep=1
   if(orig < 1) orig=nT
   #
   if(length(newxt) > 0){
      if(!is.matrix(newxt))newxt=as.matrix(newxt)
      ### calculate predictions
      h1=dim(newxt)[1]
      hstep=min(h1,hstep)
      nzt=as.matrix(zt[1:orig,])
# xt=rbind(as.matrix(xt[1:orig,]),newxt)
### changed made on 7/30/2014
      xt=rbind(xt[1:orig,,drop=FALSE],newxt)
      for (i in 1:hstep){
         tmp=Ph0
         ti=orig+i
         ## VAR part
         for (i in 1:p){
            idx=(i-1)*k
            tmp=tmp+Phi[,(idx+1):(idx+k)]%*%matrix(nzt[ti-i,],k,1)
         }
         if(m > -1){
            for (j in 0:m){
               jdx=j*dx
               tmp=tmp+beta[,(jdx+1):(jdx+dx)]%*%matrix(xt[ti-j,],dx,1)
            }
         }
         nzt=rbind(nzt,c(tmp))
      }
      ### compute standard errors of predictions
      mm=VARpsi(Phi,lag=hstep)
      Si=Sig
      se=matrix(sqrt(diag(Si)),1,k)
      if(hstep > 1){
         for (i in 2:hstep){
            idx=(i-1)*k
            wk=as.matrix(mm$psi[,(idx+1):(idx+k)])
            Si=Si+wk%*%Sig%*%t(wk)
            se1=sqrt(diag(Si))
            se=rbind(se,se1)
         }
      }
      ### Print forecasts
      cat("Prediction at origin: ",orig,"\n")
      cat("Point forecasts (starting with step 1): ","\n")
      print(round(nzt[(orig+1):(orig+hstep),],5))
      cat("Corresponding standard errors: ","\n")
      print(round(se[1:hstep,],5))
      
   }
   else{
      cat("Need new data for input variables!","\n")
   }
   #
}
############################
"VARXorder" <- function(x,exog,maxp=13,maxm=3,output=T){
   # Compute the AIC, BIC, HQ values and M-stat
   ##### This is a modified version of the old program in "VARorderE",
   ##### which uses the same number of data points.
   #####
   x1=as.matrix(x)
   exog=as.matrix(exog)
   nT=dim(x1)[1]
   k=dim(x1)[2]
   ksq=k*k
   if(maxp < 1)maxp=1
   nT1=dim(exog)[1]; m=dim(exog)[2]
   #
   if(nT1 > nT){
      cat("Adjustment made for different nobs:",c(nT,nT1), "\n")
   }
   if(nT > nT1){
      cat("Adjustment made for different nobs:",c(nT,nT1),"\n")
      nT=nT1
   }
   ###
   aic=matrix(0,maxp+1,maxm+1)
   bic=aic; hq=aic
   for (mm in 0:maxm){
      ### start with VAR(0) model, which uses just the sample means.
      isto=mm+1
      y=x1[isto:nT,]
      xm=rep(1,rep(nT-mm))
      for (i in 0:mm){
         xm=cbind(xm,exog[(isto-i):(nT-i),])
      }
      xm=as.matrix(xm)
      xpx=crossprod(xm,xm)
      xpxi=solve(xpx)
      xpy=t(xm)%*%y
      beta=xpxi%*%xpy
      y=y-xm%*%beta
      #
      s=t(y)%*%y/(nT-mm)
      enob=nT-mm
      c1=log(det(s))
      aic[1,mm+1]=c1+2*k*mm/enob
      bic[1,mm+1]=c1+log(enob)*k*m/enob
      hq[1,mm+1]=c1+2*log(log(enob))*k*m/enob
      #
      for (p in 1:maxp){
         ist=max(mm+1,p+1)
         enob=nT-ist+1
         y=as.matrix(x1[ist:nT,])
         xmtx=rep(1,enob)
         for (i in 0:mm){
            xmtx=cbind(xmtx,exog[(ist-i):(nT-i),])
          }
         for (j in 1:p){
            xmtx=cbind(xmtx,x1[(ist-j):(nT-j),])
          }
         xm1=as.matrix(xmtx)
         xpx = t(xm1)%*%xm1
         xpxinv=solve(xpx)
         xpy=t(xm1)%*%y
         beta=xpxinv%*%xpy
         resi=y-xm1%*%beta
         sse=(t(resi)%*%resi)/enob
         #print(paste("For p = ",p,"residual variance is", sse))
         d1=log(det(sse))
         npar=p*ksq+k*(mm+1)
         aic[p+1,mm+1]=d1+(2*npar)/enob
         bic[p+1,mm+1]=d1+(log(enob)*npar)/enob
         hq[p+1,mm+1]=d1+(2*log(log(enob))*npar)/enob
      }
      #end of for (mm in 0:maxm)
   }

 ind.min <- function(x){
   if(!is.matrix(x))x=as.matrix(x)
   r=dim(x)[1]
   c=dim(x)[2]
   xm=min(x)
   COntin=TRUE
   while(COntin){
      for(j in 1:c){
         idx=c(1:r)[x[,j]==xm]
         jj=j
         if(length(idx) > 0){
            ii=c(1:r)[x[,jj]==xm][1]
            kdx=c(ii,jj)
            ##print(kdx)
            COntin=FALSE
         }
      }
    }
    kdx
  }
   ## selection
   aicor=ind.min(aic)-1; bicor=ind.min(bic)-1;hqor=ind.min(hq)-1
   if(output){
      cat("selected order(p,s): aic = ",aicor,"\n")
      cat("selected order(p,s): bic = ",bicor,"\n")
      cat("selected order(p,s): hq = ",hqor,"\n")
   }
   
   VARXorder<-list(aic=aic,aicor=aicor,bic=bic,bicor=bicor,hq=hq,hqor=hqor)
 }

#### Regression model with time series errors (Multivariate case)
"REGts" <- function(zt,p,xt,include.mean=T,fixed=NULL,par=NULL,se.par=NULL,details=F){
   ## Fit a multivariate regression model with time series errors
   ### VAR model only.
   ### obtain preliminary estimation if needed.
   if(!is.matrix(zt))zt=as.matrix(zt)
   if(!is.matrix(xt))xt=as.matrix(xt)
   nT=dim(zt)[1];  k <- dim(zt)[2]
   kx=dim(xt)[2]
   if(p < 0)p=0
   if(length(par) < 1){
      m1=Mlm(zt,xt,constant=include.mean,output=F)
      par=c(m1$beta)
      r1=dim(m1$beta)[1]
      se.par=c(m1$se.beta)
      res=m1$residuals
      if(p > 0){
         m2=VAR(res,p,include.mean=F,output=F)
         par1=c(m2$coef); par=c(par,par1)
         se1=c(m2$secoef); se.par=c(se.par,se1)
         r2=dim(m2$coef)[1]
      }
      if(length(fixed) < 1){
         fixed=matrix(1,(r1+r2),k)
      }
    }
   r1=kx
   if(include.mean){r1=r1+1; kx=kx+1}
   r2=p*k
   cat("Number of parameters: ",length(par),"\n")
   cat("initial estimates: ",par,"\n")
   ### Set up lower and upper bounds
   lowerBounds=par; upperBounds=par
   npar=length(par)
   mult=1.5
   if(npar > 10)mult=1.3
   for (j in 1:npar){
      lowerBounds[j] = par[j]-mult*se.par[j]
      upperBounds[j] = par[j]+mult*se.par[j]
   }
   cat("Par. Lower-bounds: ",lowerBounds,"\n")
   cat("Par. Upper-bounds: ",upperBounds,"\n")

 RegXmtx <- function(zt,xt,p,par,include.mean,fixed){
  nT <- dim(zt)[1]; k <- dim(zt)[2]; kx <- dim(xt)[2]
   xmtx=xt; ist=p+1
   if(include.mean){
      xmtx=cbind(rep(1,nT),xmtx)
      kx=kx+1
    }
   ## setup parameter matrix for regressors
   beta=matrix(0,kx,k)
   icnt=0
   for(i in 1:k){
      idx=c(1:kx)[fixed[1:kx,i] > 0]
      if(length(idx) > 0){
         ii=length(idx)
         beta[idx,i]=par[(icnt+1):(icnt+ii)]
         icnt=icnt+ii
       }
     }
   res=zt-as.matrix(xmtx)%*%beta
   resi=res; tsxm=NULL
   if(p > 0){
      tsxm=res[(ist-1):(nT-1),]
      if(p > 1){
         for (j in 2:p){
            tsxm=cbind(tsxm,res[(ist-j):(nT-j),])
          }
       }
      #setup time-series parameter
      kp=p*k
      Phi=matrix(0,kp,k)
      fix1=fixed[(kx+1):(kx+kp),]
      for (i in 1:k){
         idx=c(1:kp)[fix1[,i] > 0]
         jj=length(idx)
         if(jj > 0){
            Phi[idx,i]=par[(icnt+1):(icnt+jj)]
            icnt=icnt+jj
         }
       }
      resi=res[ist:nT,]-as.matrix(tsxm)%*%Phi
     }
   RegXmtx <- list(xmtx=xmtx,residuals=resi)
  }
 
lRegts <- function(par,zt=zt,xt=xt,p=p,include.mean=include.mean,fixed=fixed){
   # compute the log-likelihood function of a REGts model
   nT <- dim(zt)[1];    k <- dim(zt)[2]
   nobe=nT-p
   m1=RegXmtx(zt,xt,p,par,include.mean,fixed)
   resi=m1$residuals
   ##
   ## evaluate log-likelihood function
   resi=as.matrix(resi)
   sig=t(resi)%*%resi/nobe
   ll=dmvnorm(resi,rep(0,k),sig)
   lRegts=-sum(log(ll))
  }

   # Estimate Parameters and Compute Numerically Hessian:
   if(details){
      fit = nlminb(start = par, objective = lRegts, zt=zt,xt=xt,p=p,include.mean=include.mean,fixed=fixed, 
              lower = lowerBounds, upper = upperBounds, control = list(trace=3))}
   else{
      fit = nlminb(start = par, objective = lRegts, zt=zt,xt=xt,p=p,include.mean=include.mean,fixed=fixed,
          control=list(step.min=0.2,step.max=0.5), lower = lowerBounds, upper = upperBounds)
      #fit=optim(par,lRegts,method=c("L-BFGS-B"),lower=lowerBounds,upper=upperBounds,hessian=TRUE)
   }
   epsilon = 0.0001 * fit$par
   Hessian = matrix(0, ncol = npar, nrow = npar)
   for (i in 1:npar) {
      for (j in 1:npar) {
         x1 = x2 = x3 = x4  = fit$par
         x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
         x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
         x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
         x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
         Hessian[i, j] = (lRegts(x1,zt=zt,xt=xt,p=p,include.mean=include.mean,fixed=fixed)
                         -lRegts(x2,zt=zt,xt=xt,p=p,include.mean=include.mean,fixed=fixed)
                         -lRegts(x3,zt=zt,xt=xt,p=p,include.mean=include.mean,fixed=fixed)
                         +lRegts(x4,zt=zt,xt=xt,p=p,include.mean=include.mean,fixed=fixed))/
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
   ### Prepare parameters for printing
   beta=NULL; se.beta=NULL; Phi=NULL; se.Phi=NULL; icnt=0
   ### Regressor parameters
   if(r1 > 0){
      beta=matrix(0,r1,k)
      se.beta=matrix(1,r1,k)
      for (i in 1:k){
         idx=c(1:r1)[fixed[1:r1,i] > 0]
         ii=length(idx)
         if(ii > 0){
            beta[idx,i]=est[(icnt+1):(icnt+ii)]
            se.beta[idx,i]=se.coef[(icnt+1):(icnt+ii)]
            icnt=icnt+ii
          }
       }
      cat("======= ","\n")
      cat("Coefficient matrix for constant + exogenous variable","\n")
      cat("Estimates: ","\n")
      print(round(t(beta),3))
      cat("Standard errors: ","\n")
      print(round(t(se.beta),3))
     }
   ### VAR parameters
   if(r2 > 0){
      Phi=matrix(0,r2,k)
      se.Phi=matrix(1,r2,k)
      fix1=fixed[(r1+1):(r1+r2),]
      for (i in 1:k){
         idx=c(1:r2)[fix1[1:r2,i] > 0]
         ii=length(idx)
         if(ii > 0){
            Phi[idx,i]=est[(icnt+1):(icnt+ii)]
            se.Phi[idx,i]=se.coef[(icnt+1):(icnt+ii)]
            icnt=icnt+ii
          }
       }
      cat("VAR coefficient matrices: ","\n")
      for (i in 1:p){
         kdx=(i-1)*k
         cat("AR(",i,") coefficient: ","\n")
         phi=t(Phi[(kdx+1):(kdx+k),])
         print(round(phi,3))
         cat("standard errors:","\n")
         sephi=t(se.Phi[(kdx+1):(kdx+k),])
         print(round(sephi,3))
       }
    }
   ### compute the residuals
   m1=RegXmtx(zt,xt,p,est,include.mean,fixed)
   resi=m1$residuals
   sig=t(resi)%*%resi/(nT-p)
   cat("Residual Covariance matrix: ","\n")
   print(sig,digits=4)
   d1=log(det(sig))
   aic=d1+2*npar/(nT-p)
   bic=d1+log(nT-p)*npar/(nT-p)
   cat("============","\n")
   cat("Information criteria: ","\n")
   cat("AIC: ",aic,"\n")
   cat("BIC: ",bic,"\n")
   #
   coef=rbind(beta,Phi)
   se.coef=rbind(se.beta,se.Phi)
      
   REGts <- list(data=zt,xt=xt,aror=p,include.mean=include.mean,Phi=t(Phi),se.Phi=t(se.Phi),
   beta=t(beta),se.beta=t(se.beta),residuals=resi,Sigma=sig,coef=coef,se.coef=se.coef)
}


#### Refinement
"refREGts" <- function(m1,thres=1.0){
   zt=m1$data; xt=m1$xt; p=m1$aror; include.m=m1$include.mean
   coef=m1$coef; se.coef=m1$se.coef
   k=dim(zt)[2]
   n1=dim(coef)[1]
   kx=dim(xt)[2]
   if(include.m)kx=kx+1
   icnt=0
   par=NULL; separ=NULL
   ### locate significant beta-parameters
   fix1=matrix(0,kx,k)
   for (i in 1:k){
      tt=coef[1:kx,i]/se.coef[1:kx,i]
      idx=c(1:kx)[abs(tt) > thres]
      ii=length(idx)
      if(ii > 0){
         par=c(par,coef[idx,i])
         separ=c(separ,se.coef[idx,i])
         fix1[idx,i]=1
         icnt=icnt+ii
      }
   }
   ## locate the significant VAR parameters
   r2=n1-kx
   Phi=coef[(kx+1):n1,]
   sePhi=se.coef[(kx+1):n1,]
   fix2=matrix(0,r2,k)
   for (i in 1:k){
      tt=Phi[,i]/sePhi[,i]
      jdx=c(1:r2)[abs(tt) > thres]
      jj=length(jdx)
      if(jj > 0){
         par=c(par,Phi[jdx,i])
         separ=c(separ,sePhi[jdx,i])
         icnt=icnt+jj
         fix2[jdx,i]=1
      }
   }
   fix=rbind(fix1,fix2)
   ####print(fix)
   
   mm = REGts(zt,p,xt,include.mean=include.m,fixed=fix,par=par,se.par=separ)
   coef=mm$coef;se.coef=mm$se.coef;Phi=t(mm$Phi);se.Phi=t(mm$se.Phi)
   beta=t(mm$beta);se.beta=t(mm$se.beta);resi=mm$residuals;sig=mm$Sigma
   
   refREGts <- list(data=zt,xt=xt,aror=p,include.mean=include.m,Phi=t(Phi),se.Phi=t(se.Phi),
   beta=t(beta),se.beta=t(se.beta),residuals=resi,Sigma=sig,coef=coef,se.coef=se.coef)
}


"Mlm" <- function(y,z,constant=TRUE,output=TRUE){
   # This program performs multivariate linear regression analysis.
   # z: design matrix
   # constant: switch for the constant term of the regression model
   # y: dependent variables
   ## Model is y = z%*%beta+error
   z=as.matrix(z)
   n=nrow(z)
   nx=ncol(z)
   zc=z
   if (constant) zc=cbind(rep(1,n),z)
   p=ncol(zc)
   y=as.matrix(y)
   ztz=t(zc)%*%zc
   zty=t(zc)%*%y
   ztzinv=solve(ztz)
   beta=ztzinv%*%zty
   res=y-zc%*%beta
   sig=t(res)%*%res/(n-p)
   co=kronecker(sig,ztzinv)
   sd=sqrt(diag(co))
   se.beta=matrix(sd,nrow(beta),ncol(beta))
   #
   if(output){
      print("LSE of parameters")
      print("  est   s.d.   t-ratio    prob")
      par=c(beta)
      deg=n-p
      iend=nrow(beta)*ncol(beta)
      tmp=matrix(0,iend,4)
      for (i in 1:iend){
         tt=par[i]/sd[i]
         pr=2*(1-pt(abs(tt),deg))
         tmp[i,]=c(par[i],sd[i],tt,pr)
       }
      print(tmp,digits=3)
    }
   Mlm <- list(beta=beta,se.beta=se.beta,residuals=res,sigma=sig)
 }

#### Missing value programs
####
"Vmiss" <- function(zt,piwgt,sigma,tmiss,cnst=NULL,output=T){
   ## Estimate the missing values of a vector time series.
   ## zt: T-by-k time series
   ## piwgt: pi-weight matrices of the model: k-by-N matrix
   ### [pi1,pi2,pi3, ...]. The number of columns determines the loop used in
   ### estimating the missing value.
   ## tmiss: time index of missing values
   ## sigma: k-by-k matrix of the innovations
   ##
   ### cnst: k-by-1 vector of constant.
   ###      cnst is Ph0 for VAR model and [theta(1)]^{-1}*Ph0 for VARMA models.
   ###
   if(!is.matrix(zt))zt=as.matrix(zt)
   if(!is.matrix(piwgt))piwgt=as.matrix(piwgt)
   if(!is.matrix(sigma))sigma=as.matrix(sigma)
   m1=eigen(sigma)
   va=sqrt(m1$values); P=m1$vectors; di=diag(1/va)
   Sroot=P%*%di%*%t(P)  # square root matrix of sigma
   k=dim(zt)[2]; nT=dim(zt)[1]
   if(length(cnst) < 0){
      cnst=matrix(0,1,k)}
   else{
      cnst=matrix(cnst,1,k)
   }
   #
   lpi=dim(piwgt)[2]
   lags=floor(lpi/k)
   ## setup the multivariate linear regression
   Y=NULL; X=NULL; Tpiwgt=t(piwgt)
   zt1=zt; zt1[tmiss,]=rep(0,k)
   iend=min(nT,(tmiss+lags))
   if((tmiss > 1) && (tmiss < nT)){
      # estimate the missing values
      wk=matrix(0,1,lpi)
      icnt=0
      jend=min(tmiss-1,lags)
      for (j in 1:jend){
         wk[1,(icnt+1):(icnt+k)]=zt1[tmiss-j,]
         icnt=icnt+k
      }
      yt=wk%*%Tpiwgt+cnst
      xt=diag(rep(1,k))
      Y=Sroot%*%t(yt)
      X=Sroot%*%xt
      #
      Tmax=min(nT,tmiss+lags)
      iend=min(lags,Tmax-tmiss)
      for (i in 1:iend){
         wk[1,]=c(zt1[tmiss+i-1,],wk[1,1:(lpi-k)])
         yt=zt1[(tmiss+i),]-wk%*%Tpiwgt-cnst
         kst=(i-1)*k
         xt=piwgt[,(kst+1):(kst+k)]
         Y=rbind(Y,Sroot%*%t(yt))
         X=rbind(X,Sroot%*%xt)
      }
      xpx=t(X)%*%X
      xpy=t(X)%*%Y
      xpxi=solve(xpx)
      est=xpxi%*%xpy
      if(output){
         cat("Estimate of missing value at time index",tmiss,"\n")
         print(round(est,5))
       }
   }
   return(est)
 }

"Vpmiss" <- function(zt,piwgt,sigma,tmiss,mdx,cnst=NULL,output=T){
   ## Estimate a PARTIALLY missing observation of a vector time series.
   ## See Vmiss for variable descriptions
   ### mdx: a k-dimensional vector to locating missing components.
   ### mdx[i] = 0 if the i-th component is missing;
   ### mdx[i] = 1 if the -th component is observed.
   ###
   if(!is.matrix(zt))zt=as.matrix(zt)
   if(!is.matrix(piwgt))piwgt=as.matrix(piwgt)
   if(!is.matrix(sigma))sigma=as.matrix(sigma)
   k=dim(zt)[2]; nT=dim(zt)[1]
   if(length(mdx) < 1)mdx=rep(0,k)
   miss=c(1:k)[mdx==0]
   nmiss=c(1:k)[mdx==1]
   nm=length(miss)
   arrange=c(miss,nmiss)
   ARR=diag(rep(1,k))
   ARR=ARR[,arrange]
   ARRi=solve(ARR)
   ### re-arrange the Sigma matrix
   sig=ARR%*%sigma%*%t(ARR)
   m1=eigen(sig)
   va=sqrt(m1$values); P=m1$vectors; di=diag(1/va)
   Sroot=P%*%di%*%t(P)  # square root matrix of sigma
   #
   if(length(cnst) < 0){
      cnst=matrix(0,1,k)}
   else{
      cnst=matrix(cnst,1,k)
   }
   #
   ncnst=t(ARR%*%t(cnst))
   lpi=dim(piwgt)[2]
   lags=floor(lpi/k)
   ## compute the transformed pi-weight matrices
   npiwgt=piwgt
   for (i in 1:lags){
      icnt=(i-1)*k
      tmp=piwgt[,(icnt+1):(icnt+k)]
      tmp1=ARR%*%tmp%*%ARRi
      npiwgt[,(icnt+1):(icnt+k)]=tmp1
   }
   ##### Partially missing can occur at any time point
   ##### This part is different from "Vmiss.R". However, estimation of
   ##### missing values at the beginning of a time series may not be efficient.
   ## setup the multivariate linear regression
   Y=NULL; X=NULL; Tpiwgt=t(npiwgt)
   zt1=zt[,arrange]; zt1[tmiss,1:nm]=rep(0,nm)
   yobs=rep(0,k)
   if(k > nm)yobs=matrix(zt1[tmiss,(nm+1):k],(k-nm),1)
   zt1[tmiss,]=rep(0,k)
   iend=min(nT,(tmiss+lags))
   # estimate the missing values
   wk=matrix(0,1,lpi)
   icnt=0
   jend=min(tmiss-1,lags)
   if(jend > 0){
      for (j in 1:jend){
         wk[1,(icnt+1):(icnt+k)]=zt1[tmiss-j,]
         icnt=icnt+k
       }
    }
   yt=wk%*%Tpiwgt+ncnst
   xt=diag(rep(1,k))
   Y1=Sroot%*%t(yt)
   X1=Sroot%*%xt
   if(k > nm){
      Y=Y1-X1[,(nm+1):k]%*%yobs
     }
   else{
      Y=Y1
    }
   X=X1[,1:nm]
   if(nm == 1)X=matrix(X1[,1],k,1)
   Tmax=min(nT,tmiss+lags)
   iend=min(lags,Tmax-tmiss)
   if(iend > 0){
      for (i in 1:iend){
         if(lpi > k){
           wk[1,]=c(zt1[tmiss+i-1,],wk[1,1:(lpi-k)])
           }
           else{
            wk[1,]=zt1[tmiss+i-1,]
            }
         yt=zt1[(tmiss+i),]-wk%*%Tpiwgt-ncnst
         kst=(i-1)*k
         xt=npiwgt[,(kst+1):(kst+k)]
         Y1=Sroot%*%t(yt)
         X1=Sroot%*%as.matrix(xt)
         if(k > nm){
            Y=rbind(Y,Y1-X1[,(nm+1):k]%*%yobs)
          }
         else{
            Y=rbind(Y,Y1)
          }
         if(nm==1){
            X1=matrix(X1[,1],k,1)
            X=rbind(X,X1)
          }
         else{
            X=rbind(X,X1[,1:nm])}
      }
   }
   xpx=t(X)%*%X
   xpy=t(X)%*%Y
   xpxi=solve(xpx)
   est=xpxi%*%xpy
   if(output){
      cat("Estimate of missing value at time index",tmiss,"\n")
      cat("The missing idicator is: ",mdx,"\n")
      print(round(est,5))
    }
   return(est)
 }

"SCMid" <- function(zt,maxp=5,maxq=5,h=0,crit=0.05,output=FALSE){
   ### Identify SCMs for a given k-dimensional time series zt.
   if(!is.matrix(zt))zt=as.matrix(zt)
   nT=dim(zt)[1]
   k=dim(zt)[2]
   nar=maxp+1; nma=maxq+1
   zeroTbl=matrix(0,nar,nma)
   diagDif=matrix(0,nar,nma)
   for (m in 0:maxp){
      for (j in 0:maxq){
         ##h=maxq-1-j
         ist=m+j+h+2
         ### setup the Ymt matrix
         Ymt=zt[ist:nT,]
         if(m > 0){
            for (i in 1:m){
               Ymt=cbind(Ymt,zt[(ist-i):(nT-i),])
            }
         }
         Ymt=as.matrix(Ymt)
         k1=dim(Ymt)[2]
         ### setup the Yht-matrix and denote it by Pt.
         Pt=zt[(ist-j-1):(nT-j-1),]
         if(m > 0){
            for (i in 1:m){
               Pt=cbind(Pt,zt[(ist-1-i-j):(nT-1-i-j),])
            }
         }
         if(h > 0){
            for (i in 1:h){
               Pt=cbind(Pt,zt[(ist-1-j-m-i):(nT-1-j-m-i),])
            }
         }
         Pt=as.matrix(Pt)
         ##print(c(dim(Ymt),dim(Pt)))
         m1=cancor(Ymt,Pt)
         corsq=m1$cor^2
         #
         if(output){
            cat("For (m,j) = ",c(m,j),"\n")
            cat("Squares of canonical correlations: ","\n")
            print(round(sort(corsq),3))
         }
         ### compute the variance of canonical correlations
         dsq=rep(1,k1)
         if(j > 0){
            xM=as.matrix(m1$xcoef)
            yM=as.matrix(m1$ycoef)
            ## Normalization
            for (kk in 1:k1){
               xM[,kk]=xM[,kk]/sqrt(sum(xM[,kk]^2))
               yM[,kk]=yM[,kk]/sqrt(sum(yM[,kk]^2))
            }
            xT=Ymt%*%xM
            yT=Pt%*%yM[,1:k1]
            for (jj in 1:k1){
               d1=1.0
               m1a=acf(xT[,jj],lag.max=j,plot=F)
               m1b=acf(yT[,jj],lag.max=j,plot=F)
               for (ij in 2:(j+1)){
                  d1=d1+2*m1a$acf[ij]*m1b$acf[ij]
               }
               dsq[jj]=d1
               ## end of jj-loop
            }
            ## end of if(j > 0)
         }
         #### Since dsq is only valid for stationary series, we set an upper limit.
         chk=qnorm(1-crit/2)
         chk1=chk^2/nT
         idx=c(1:k1)[corsq > chk1 ]
         dsq[idx]=1
         ###
         if(output){
            cat("Variance of canonical correlations","\n")
            print(round(dsq,3))
            cat("Test results: ","\n")
         }
         ### Perform tests to check the number of SCMs at the (m,j)-position
         icnt=0
         tst=0
         k2=(j+1)*k
         k2=min(k1,k2)
         for (kk in 1:k2){
            tmp=corsq[k1-kk+1]/dsq[k1-kk+1]
            if(tmp >= 1)tmp=0.999
            tmp=log(1-tmp)
            tst=tst-(nT-m-j)*tmp
            deg=kk*(kk+h*k)
            pv=1-pchisq(tst,deg)
            if(output)print(c(tst,deg,pv),digits=3)
            if(pv >= crit)icnt=icnt+1
         }
         zeroTbl[(m+1),(j+1)]=icnt
       }
     }
   #### Print the output tables
   cat("Column: MA order","\n")
   cat("Row   : AR order","\n")
   cat("Number of zero canonical correlations","\n")
   colnames(zeroTbl) <- c(c(0:maxq))
   rownames(zeroTbl) <- c(c(0:maxp))
   printCoefmat(zeroTbl)
   diagDif=zeroTbl
   for (i in 1:maxp){
      for (j in 1:maxq){
         diagDif[(i+1),(j+1)]=min(zeroTbl[(i+1),(j+1)]-zeroTbl[i,j],k)
       }
    }
   cat("Diagonal Differences: ","\n")
   colnames(diagDif) <- c(c(0:maxq))
   rownames(diagDif) <- c(c(0:maxp))
   printCoefmat(diagDif)
   
   SCMid <- list(Nmtx=zeroTbl,DDmtx=diagDif)
 }

##### Second-stage of specifiction ##################################################
"SCMid2" <- function(zt,maxp=2,maxq=2,h=0,crit=0.05,sseq=NULL){
### Identify details of specified SCMs. This is a second-step specification.
### sseq denotes the sequence of orders (m,j) for searching SCMs.
####
if(!is.matrix(zt))zt=as.matrix(zt)
nT=dim(zt)[1]
k=dim(zt)[2]
nar=maxp+1; nma=maxq+1; Order=matrix(0,k,2)
cmax=max((maxp+1)*k*k*k,100)
maxpq=maxp+maxq
if(length(sseq) < 1){
## set default sequence to find SCMs
 sseq=c(0,0)
 for (ell in 1:maxpq){
   for (jj in 0:ell){
     m=ell-jj; j=jj
     if( (m <= maxp) && (j <= maxq)){sseq=rbind(sseq,c(m,j))}
     }
    }
   }
tcases=dim(sseq)[1]  ## The total number of SCM orders to be tested.
## wkspace: storage space to store genuine eigenvectors. 
wkspace=matrix(0,tcases,cmax)
Nscm=rep(0,tcases)
## Tmx: Transformation matrix.
Tmx=NULL
## Jcnt: counts of linearly independent SCMs
Jcnt=0
##
for (nc in 1:tcases){
if(Jcnt >= k)break
 m=sseq[nc,1];j=sseq[nc,2]
 ist=m+j+h+2
### setup the Ymt matrix
 Ymt=zt[ist:nT,]
 if(m > 0){
  for (i in 1:m){
   Ymt=cbind(Ymt,zt[(ist-i):(nT-i),])
   }
  }
Ymt=as.matrix(Ymt)
k1=dim(Ymt)[2]
### setup the Yht-matrix and denote it by Pt.
Pt=zt[(ist-j-1):(nT-j-1),]
 if(m > 0){
  for (i in 1:m){
   Pt=cbind(Pt,zt[(ist-1-i-j):(nT-1-i-j),])
   }
  }
 if(h > 0){
  for (i in 1:h){
   Pt=cbind(Pt,zt[(ist-1-j-m-i):(nT-1-j-m-i),])
   }
  }
Pt=as.matrix(Pt)
mcan=cancor(Ymt,Pt)
corsq=mcan$cor^2
#
cat("For (pi,qi) = (",m,",",j,")","\n")

### compute the variance of canonical correlations
dsq=rep(1,k1)
xM=as.matrix(mcan$xcoef)
yM=as.matrix(mcan$ycoef)
## Normalization
for (kk in 1:k1){
  xM[,kk]=xM[,kk]/sqrt(sum(xM[,kk]^2))
  yM[,kk]=yM[,kk]/sqrt(sum(yM[,kk]^2))
  }
xT=Ymt%*%xM
yT=Pt%*%yM[,1:k1]
for (jj in 1:k1){
  d1=1.0
  if(j > 0){
   m1a=acf(xT[,jj],lag.max=j,plot=F)
   m1b=acf(yT[,jj],lag.max=j,plot=F)
   for (ij in 2:(j+1)){
    d1=d1+2*m1a$acf[ij]*m1b$acf[ij]
   }
  }
  dsq[jj]=d1
 }

#### Since dsq is only valid for stationary series, we set an upper limit.
chk=qnorm(1-crit/2)
chk1=chk^2/nT
idx=c(1:k1)[corsq > chk1 ]
dsq[idx]=1
###
cat("Tests:","\n")
re=NULL
### Perform tests to check the number of SCMs at the (m,j)-position
icnt=0
tst=0
k2=(j+1)*k
k2=min(k1,k2)
for (kk in 1:k2){
ik=k1-kk+1
tmp=corsq[ik]/dsq[ik]
if(tmp >= 1)tmp=corsq[ik]
tmp=log(1-tmp)
tst=tst-(nT-m-j)*tmp
deg=kk*(kk+h*k)
pv=1-pchisq(tst,deg)
re=rbind(re,c(corsq[ik],dsq[ik],tst,deg,pv))
if(pv >= crit)icnt=icnt+1
}
re1=round(re,3)
colnames(re1) <- c("Eigvalue","St.dev","Test","deg","p-value")
print(re1)
#
cat("Summary:","\n")
cat("Number of SCMs detected: ",icnt,"\n")
Nscm[nc]=icnt; wk1 = NULL; n12=0
#####cat("Nscm[nc]", Nscm[nc],"\n")
if(icnt > 0){
 ldx=rev(c((k1-icnt+1):k1))
 wk1=as.matrix(xM[,ldx])
###
if((nc > 1) && (Jcnt > 0)){
  ## checking for newly detected SCM, if any.
  wk2=NULL  ## contains all the implied SCMs
   ### obtain the SCMs found before.
    for (j1 in 1:(nc-1)){
     wk3=NULL
     ndup=min(m-sseq[j1,1],j-sseq[j1,2])
     m1=sseq[j1,1]; wcnt=Nscm[j1]
     if((ndup >= 0) && (wcnt > 0)){
      leng=(m1+1)*k*wcnt
      wvector=matrix(wkspace[j1,1:leng],(m1+1)*k,wcnt)
      cat("wvector","\n")
      print(wvector,digits=3)
      if(ndup == 0){
       if(m1==m){wk3=cbind(wk3,wvector)}
        else{tmp=rbind(wvector,matrix(0,(m-m1)*k,wcnt))
             wk3=cbind(wk3,tmp)
             }
       }
      if(ndup > 0){
       tmp=rbind(wvector,matrix(0,(m-m1)*k,wcnt))
       wk3=cbind(wk3,tmp)
       jdim=dim(wk3)[1]; tmp1=tmp
       for (ell in 1:ndup){
         tmp1=rbind(matrix(0,k,wcnt),tmp1[1:(jdim-k),])
         wk3=cbind(wk3,tmp1)
        }
       }
      }
     wk2=cbind(wk2,wk3)
     }
     n12=dim(wk2)[2]
     if(icnt <= n12){
      cat("No new SCM found.","\n")}
      else{
       chknew=cancor(wk1,wk2)
       xcoef=as.matrix(chknew$xcoef)
       dd=t(xcoef)%*%xcoef
       dd=sqrt(diag(dd))
       for (ii in 1:icnt){
        xcoef[,ii]=xcoef[,ii]/dd[ii]
        }
      wk=wk1%*%xcoef[,(n12+1):icnt]
      ngenu=icnt-n12
      cat("The number of newly found SCMs: ",ngenu,"\n")
      wk1=wk
      cat("Vectors: ","\n")
      print(wk1,digits=3)
      }
##  Checking for possible exchangeable SCM 
   if(ngenu > 0){
    iexch = 0       
    for (j1 in 1:(nc-1)){
     m1=sseq[j1,1]; j2=sseq[j1,2]
     if((m1+j2)==(m+j)){
       wcnt=Nscm[j1]; leng=(m1+1)*k*wcnt
       wvector = matrix(wkspace[j1,1:leng],(m1+1)*k,wcnt)
       if((wcnt > 0) && (j2 > j)){
           wk2 = rbind(wvector,matrix(0,(m-m1)*k,wcnt)); wk3=wk1}
        else{ij=dim(wk1)[1]; wk2 = wvector[1:ij,]; wk3=wk1}
       mchg=cancor(wk3,wk2)
       ich=0
       h11=min(wcnt,icnt)
       for (ij in 1:h11){
        if(mchg$cor[ij] > 0.8){ich=ich+1}
        }
       if(ich > 0){
        cat("Exchangeable SCM found with order: ", c(m1,j2),"\n")
         iexch=iexch+ich
         }
        } # end of (m1+j2) == m+j
      }  ##end of for-loop
     } ## end of if(ngenu > 0)   
   }### end of (if (nc > 1) && (Jcnt > 0))
    newSCM=icnt-n12
    if(newSCM > 0){
      for (ii in 1:newSCM){
      Order[Jcnt+ii,1]=m; Order[Jcnt+ii,2]=j
      }
     }
  } ## end of if(icnt > 0)
    if(n12 > 0){icnt=icnt-n12}
##### Update the new information found, if any
    Nscm[nc]=icnt
    Jcnt=Jcnt+icnt
    leng=(m+1)*k*icnt
    if(leng > 0){wkspace[nc,1:leng] = c(wk1)}
    Tmx=cbind(Tmx,wk1[1:k,])
  } ### end of for-loop
   cat("SUMMARY:","\n")
   cat("Overall model: ",apply(Order,2,max),"\n")
   cat("Orders of SCM: ","\n")
   print(Order)
    cat("Transformation Matrix: ","\n")
    print(Tmx,digits=3)

  SCMid2 <- list(Tmatrix = t(Tmx),SCMorder=Order)
 }

"SCMmod" <- function(order,Ivor,output){
   ## Parameter specification for a given set of SCMs and indicator of T-matrix.
   ## order is a k-by-2 matrix of the orders of SCMs.
   ## Ivor is a k-dimensional vector with position of 1.
   if(!is.matrix(order))order=as.matrix(order)
   k=dim(order)[1]
   p=max(order[,1])
   q=max(order[,2])
   kp=k*p; kq=k*q
   if(output) cat("VARMA order (p,q) = (",p,",",q,")","\n")
   phi=NULL; theta=NULL
   if(p > 0)phi=matrix(0,k,kp)
   if(q > 0)theta=matrix(0,k,kq)
   ### 2: estimation; 0: fixed at zero; 1: fixed at 1.
   for (i in 1:k){
      pi=order[i,1]; qi=order[i,2]
      if(pi > 0){
         for (j in 1:pi){
            jst=k*(j-1)
            phi[i,(jst+1):(jst+k)]=2
         }
      }
      if(qi > 0){
         for (j in 1:qi){
            jst=k*(j-1)
            theta[i,(jst+1):(jst+k)]=2
         }
      }
   }
   ## check for redundant parameters
   for (i in 1:k){
      pi=order[i,1];qi=order[i,2]
      for (j in 1:i){
         pj=order[j,1];qj=order[j,2]
         mm=min(pi-pj,qi-qj)
         if(mm > 0){
            for (jj in 1:mm){
               jdx=(jj-1)*k+j
               theta[i,jdx]=0
            }
         }
         ## finish redundant parameters
      }
      # finish the AR and MA parameter sepcification
   }
   if(output){
      cat("VAR matrices: ","\n")
      print(phi)
      cat("VMA matrices: ","\n")
      print(theta)
   }
   #
   Tmtx=matrix(2,k,k)
   ### finding the pivotal positions
   for (i  in 1:k){
      k1=Ivor[i]
      Tmtx[i,k1]=1
   }
   # row-by-row
   for (i in 2:k){
      pi=order[i,1]; qi=order[i,2]
      for (j in 1:(i-1)){
         k1=Ivor[j]
         pj=order[j,1];qj=order[j,2]
         mink=min(pi-pj,qi-qj)
         if(mink > -1)Tmtx[i,k1]=0
       }
    }
   if(output){
      cat("Transformation matrix: ","\n")
      print(Tmtx)
    }
   SCMmod <- list(Tmtx=Tmtx,ARpar=phi,MApar=theta)
 }

#### SCM estimation
"SCMfit" <- function(da,scms,Tdx,include.mean=T,fixed=NULL,prelim=F,details=F,thres=1.0,ref=0,SCMpar=NULL,seSCMpar=NULL){
   # Estimation of a vector ARMA model using conditional MLE (Gaussian dist)
   #  The model is specified via SCMs.
   #
   # The program is modified from Kronfit.R, July 2012. 
   #
   # When prelim=TRUE, fixed is assigned based on the results of AR approximation.
   # Here "thres" is used either prelim = TRUE or in refined estimation.
   # ref = 0 denotes not a refined estimation.
   ##
   # use mvtnorm package for multivariate normal density.
   if(!is.matrix(da))da=as.matrix(da)
   nT=dim(da)[1]; k=dim(da)[2]
   p=max(scms[,1]); q=max(scms[,2]); pq=max(p,q); kp=k*p; kq=k*q
   pq=max(p,q)
   cat("Maximum VARMA order: (",p,",",q,")","\n")
      # Obtain the parameter specifications.
      mm1=SCMmod(scms,Tdx,FALSE)
      # Step 1: assign the data and locations of parameters globally
      locTmtx <- mm1$Tmtx; locAR <- mm1$ARpar; locMA <- mm1$MApar
      #
   if(ref < 1){
      cat("Locations of estimable parameters: Transformation Matrix","\n")
      print(locTmtx)
      cat("AR parameters","\n")
      print(locAR)
      cat("MA parameters","\n")
      print(locMA)
    }
 iniSCM <- function(da,at,scms,Tdx,locTmtx,locAR,locMA,inc.mean){
   ### The output parameters are in a vector:
   #### It contains parameters equation-by-equation, including the constant term, if any.
   ##### Format:
   #### z(t) = xi0 z(t) + SUM[xii *z(t-i)] + SUM[Omegai*a(t-i)] + a(t).
   if(!is.matrix(da))da=as.matrix(da)
   if(!is.matrix(at))at=as.matrix(at)
   nT=dim(da)[1]; k=dim(da)[2]; p = 0; q=0
   ## obtain the maximum index value.
   if(length(locAR) > 0)p=floor(dim(locAR)[2]/k)
   if(length(locMA) > 0)q=floor(dim(locMA)[2]/k)
   ##   cat("iniSCM p and q: ",c(p,q),"\n")
   pq=p+q
   if(pq < 1){
      cat("The series is white noise. No estimation is needed","\n")
   }
   ## est: stores the estimates (equation 1, equation 2, etc.)
   est=NULL
   estse=NULL
   ist=pq+1
   ### Perform estimation equation-by-equation.
   for (i in 1:k){
      ## setup the variables of the contemporaneous lag if any.
      X=NULL
      loki=Tdx[i]
      Y=da[ist:nT,loki]
      ##
      if(inc.mean)X=rep(1,(nT-pq))
      ### setup contemporaneous parameters
      for(j in 1:k){
         if(j != loki){
            if(locTmtx[i,j] > 1)X=cbind(X,-da[ist:nT,j])
         }
      }
      ### setup the lagged AR variables
      pi=scms[i,1]; qi=scms[i,2]
      if(pi > 0){
         for(lag in 1:pi){
            jst=(lag-1)*k
            for (j in 1:k){
               if(locAR[i,jst+j] > 1){
                  tmp=da[(ist-lag):(nT-lag),j]
                  X=cbind(X,tmp)
               }
            }
         }
      }
      ### Next come the MA part
      ### setup the lagged MA variables
      if(qi > 0){
         for(lag in 1:qi){
            jst=(lag-1)*k
            for (j in 1:k){
               if(locMA[i,jst+j] > 1){
                  tmp=at[(ist-lag):(nT-lag),j]
                  X=cbind(X,tmp)
               }
            }
         }
       }
      ### Perform estimation
      XPX=crossprod(X,X)/nT
      XPXinv=solve(XPX)
      XPY=crossprod(X,Y)/nT
      beta=XPXinv%*%XPY
      l1=dim(XPX)[1]
      resi=Y-X%*%matrix(beta,l1,1)
      evar=crossprod(resi,resi)/(nT-p)
      est=c(est,beta)
      estse=c(estse,sqrt(diag(XPXinv)*evar/nT))
   }
   iniSCM <- list(par=est,se=estse)
 }
   # The next statement is designed for model refinement.
   if(ref < 1){
       ### Use VAR approximation to obtain approximate innovations
      m1=VARorder(da,pq+9,output=FALSE)
      porder=m1$aicor
      if(porder < 1)porder=1
      m2=VAR(da,porder,output=FALSE)
      y=da[(porder+1):nT,]
      x=m2$residuals
      m3=iniSCM(y,x,scms,Tdx,locTmtx,locAR,locMA,include.mean)
      ### SCMpar is the vector of ALL estimable parameters. 
      ####               [Some of which maybe fixed to zero.]
      SCMpar <- m3$par; seSCMpar <- m3$se 
      ### SCMpar is a vector; which stores parameters equation-by-equation.
      ##
      nr=length(SCMpar)
      ### Preliminary simplification
      if(prelim){
         fixed = rep(0,nr)
         for (j in 1:nr){
            tt=SCMpar[j]/seSCMpar[j]
            if(abs(tt) >= thres){
               fixed[j]=1
            }
            else{
               SCMpar[j]=0
            }
            # end of j-loop
         }
         # end of if(prelim)
      }
      #### fixed = 1 means "estimation"; = 0 means fixed to zero.
      if(length(fixed) < 1)fixed=rep(1,nr)
    }
   #### JJdx is kept for identification purpose.
   nr=length(SCMpar)
   JJdx=c(1:nr)[fixed==1]
   par=SCMpar[JJdx]
   separ= seSCMpar[JJdx]
   #########
   cat("Number of parameters: ",length(par),"\n")
   cat("initial estimates: ",round(par,4),"\n")
   ### Set up lower and upper bounds
   lowerBounds=par; upperBounds=par
   for (j in 1:length(par)){
      lowerBounds[j] = par[j]-2*separ[j]
      upperBounds[j] = par[j]+2*separ[j]
   }
   cat("Upper-bound: ",round(upperBounds,4),"\n")
   cat("Lower-bound: ",round(lowerBounds,4),"\n")

 LLSCM <- function(par,zt=da,scms=scms,Tdx=Tdx,SCMpar=SCMpar,JJdx=JJdx,include.mean=include.mean,fixed=fixed,locTmtx=locTmtx,locAR=locAR,locMA=locMA){
   k <- dim(zt)[2];  nT <- dim(zt)[1]
   p=max(scms[,1]); q=max(scms[,2]); pq=max(p,q); kp=k*p; kq=k*q
   Tdx <- Tdx
   SCMpar[JJdx] <- par
   ###  Assign parameters to their proper locations in the program.
   Cnt=rep(0,k)
   ### separate lag-0 coefficient matrix.
   Ph0=locTmtx; PH=NULL; TH=NULL
   if(p > 0)PH=matrix(0,k,kp) 
   if(q > 0)TH=matrix(0,k,kq)
   icnt=0
   for (i in 1:k){
      idx=NULL; jdx=NULL; kdx=NULL
      if(p > 0)idx=c(1:kp)[locAR[i,] > 1]
      if(q > 0)jdx=c(1:kq)[locMA[i,] > 1]
      # kdx denotes the number of non-zero elements in lag-0.
      kdx=c(1:k)[locTmtx[i,] > 1]
      iend=length(idx); jend=length(jdx); kend=length(kdx)
      #### icnt: parameter count
      if(include.mean){
         icnt=icnt+1
         Cnt[i]=SCMpar[icnt]
       }
      if(kend > 0){
         Ph0[i,kdx]=SCMpar[(icnt+1):(icnt+kend)]
         icnt=icnt+kend
       }
      if(iend > 0){
         PH[i,idx]=SCMpar[(icnt+1):(icnt+iend)]
         icnt=icnt+iend
       }
      if(jend > 0){
         TH[i,jdx]=SCMpar[(icnt+1):(icnt+jend)]
         icnt=icnt+jend
       }
    }
   ##### Compute the residuals 
   ###### Compute the AR and MA coefficient matrix
   Ph0i=solve(Ph0); ARc=NULL; MAc=NULL
   if(p > 0)ARc=Ph0i%*%PH
   if(q > 0)MAc=Ph0i%*%TH
   Cntc=Ph0i%*%as.matrix(Cnt,k,1)
   ## 
   ist=pq+1
   #### consider the case t from 1 to pq+1
   at=matrix((zt[1,]-Cntc),1,k)
   if(pq > 1){
      for (t in 2:pq){
         tmp=matrix((zt[t,]-Cntc),1,k)
         if(p > 0){
            for (j in 1:p){
               if((t-j) > 0){
                  jdx=(j-1)*k
                  tmp1=matrix(zt[(t-j),],1,k)%*%t(as.matrix(ARc[,(jdx+1):(jdx+k)]))
                  tmp=tmp-tmp1
               }
             }
          }
         if(q > 0){
            for (j in 1:q){
               jdx=(j-1)*k
               if((t-j)>0){
                  tmp2=matrix(at[(t-j),],1,k)%*%t(as.matrix(MAc[,(jdx+1):(jdx+k)]))
                  tmp=tmp-tmp2
               }
            }
          }
         at=rbind(at,tmp)
      }
   }
   ### for t from ist on
   ist=pq+1
   Pcnt = NULL; beta=NULL
   if(q < 1)MAc=NULL; if(p < 1)ARc=NULL
   if(include.mean)beta=matrix(Cntc,1,k)
   if(length(ARc) > 0)beta=rbind(beta,t(ARc))
   if(length(MAc) > 0)beta=rbind(beta,t(MAc))
   #
   idim=k*(p+q)
   if(include.mean){
      Pcnt=c(1)
      idim=idim+1
   }
   #
   for (t in ist:nT){
      Past=NULL
      if(p > 0){
         for (j in 1:p){
            Past=c(Past,zt[(t-j),])
         }
      }
      if(q > 0){
         for (j in 1:q){
            Past=c(Past,at[(t-j),])
         }
      }
      tmp = matrix(c(Pcnt,Past),1,idim)%*%beta
      tmp3=zt[t,]-tmp
      at=rbind(at,tmp3)
   }
   at=as.matrix(at[(ist:nT),])
   sig=crossprod(at,at)/(nT-pq)
   ll=dmvnorm(at,rep(0,k),sig)
   LLSCM=-sum(log(ll))
   LLSCM
 }

  #  Estimate Parameters and Compute Numerically Hessian:
   if(details){ 
      fit = nlminb(start = par, objective = LLSCM,zt=da,scms=scms,Tdx=Tdx,SCMpar=SCMpar,JJdx=JJdx,include.mean=include.mean,fixed=fixed,
      locTmtx=locTmtx,locAR=locAR,locMA=locMA,lower = lowerBounds, upper = upperBounds, control = list(trace=3))
   }
   else {
      fit = nlminb(start = par, objective = LLSCM,zt=da,scms=scms,Tdx=Tdx,SCMpar=SCMpar,JJdx=JJdx,include.mean=include.mean,fixed=fixed, 
       locTmtx=locTmtx,locAR=locAR,locMA=locMA,lower = lowerBounds, upper = upperBounds)
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
         Hessian[i, j] = (LLSCM(x1,zt=da,scms=scms,Tdx=Tdx,SCMpar=SCMpar,JJdx=JJdx,include.mean=include.mean,fixed=fixed,locTmtx=locTmtx,locAR=locAR,locMA=locMA)
                         -LLSCM(x2,zt=da,scms=scms,Tdx=Tdx,SCMpar=SCMpar,JJdx=JJdx,include.mean=include.mean,fixed=fixed,locTmtx=locTmtx,locAR=locAR,locMA=locMA)
                         -LLSCM(x3,zt=da,scms=scms,Tdx=Tdx,SCMpar=SCMpar,JJdx=JJdx,include.mean=include.mean,fixed=fixed,locTmtx=locTmtx,locAR=locAR,locMA=locMA)
                         +LLSCM(x4,zt=da,scms=scms,Tdx=Tdx,SCMpar=SCMpar,JJdx=JJdx,include.mean=include.mean,fixed=fixed,locTmtx=locTmtx,locAR=locAR,locMA=locMA))/
         (4*epsilon[i]*epsilon[j])
      }
   }
   # Step 6: Create and Print Summary Report:
   d1=det(Hessian)
   if(d1 < 1.0e-13){
      se.coef=rep(1,npar)
   }
   else{
      se.coef = sqrt(diag(solve(Hessian)))
   }
   tval = fit$par/se.coef
   matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
   dimnames(matcoef) = list(names(tval), c(" Estimate",
   " Std. Error", " t value", "Pr(>|t|)"))
   cat("\nCoefficient(s):\n")
   printCoefmat(matcoef, digits = 4, signif.stars = TRUE)
   SCMpar[JJdx]=fit$par
   seSCMpar[JJdx]=se.coef
   ##cat("SCMpar: ",round(SCMpar,3),"\n")
   ######### Use locTmtx, locAR, and locMA and include.mean to identify the parameters.
   # Restore estimates to the format of unconstrained case for printing.
   Cnt=rep(0,k); seCnt=rep(0,k)
   ### separate lag-0 coefficient matrix.
   Ph0=locTmtx; sePh0=matrix(0,k,k); PH=NULL; TH=NULL
   if(p > 0){
      PH=matrix(0,k,kp); sePH=matrix(0,k,kp)
   }
   else {
      PH = NULL; sePH=NULL
   }
   if(q > 0){
      TH=matrix(0,k,kq); seTH=matrix(0,k,kq)
   }
   else{
      TH=NULL; seTH=NULL
   }
   ###
   icnt=0
   for (i in 1:k){
      idx=NULL; jdx=NULL; kdx=NULL
      if(p > 0)idx=c(1:kp)[locAR[i,] > 1]
      if(q > 0)jdx=c(1:kq)[locMA[i,] > 1]
      kdx=c(1:k)[locTmtx[i,] > 1]
      iend=length(idx); jend=length(jdx); kend=length(kdx)
      #### icnt: parameter count
      if(include.mean){
         icnt=icnt+1
         Cnt[i]=SCMpar[icnt]
         seCnt[i]=seSCMpar[icnt]
      }
      if(kend > 0){
         Ph0[i,kdx]=SCMpar[(icnt+1):(icnt+kend)]
         sePh0[i,kdx]=seSCMpar[(icnt+1):(icnt+kend)]
         icnt=icnt+kend
      }
      if(iend > 0){
         ##cat("idx: ",idx,"\n")
         PH[i,idx]=SCMpar[(icnt+1):(icnt+iend)]
         sePH[i,idx]=seSCMpar[(icnt+1):(icnt+iend)]
         icnt=icnt+iend
      }
      if(jend > 0){
         TH[i,jdx]=SCMpar[(icnt+1):(icnt+jend)]
         seTH[i,jdx]=seSCMpar[(icnt+1):(icnt+jend)]
         icnt=icnt+jend
      }
      ### end of the i-loop
   }
   #########
   cat("---","\n")
   cat("Estimates in matrix form:","\n")
   if(include.mean){
      cat("Constant term: ","\n")
      cat("Estimates: ",round(Cnt,3),"\n")
   }
   cat("AR and MA lag-0 coefficient matrix","\n")
   print(round(Ph0,3))
   jcnt=0
   if(p > 0){
      cat("AR coefficient matrix","\n")
      for (i in 1:p){
         cat("AR(",i,")-matrix","\n")
         ph=PH[,(jcnt+1):(jcnt+k)]
         print(round(ph,3))
         jcnt=jcnt+k
      }
   }
   if(q > 0){
      cat("MA coefficient matrix","\n")
      icnt=0
      for (i in 1:q){
         cat("MA(",i,")-matrix","\n")
         theta=-TH[,(icnt+1):(icnt+k)]
         print(round(theta,3))
         icnt=icnt+k
      }
   }
   ##### Compute the residuals 
   ###### Compute the AR and MA coefficient matrix
   Ph0i=solve(Ph0); ARc = NULL; MAc=NULL
   if(p > 0)ARc=Ph0i%*%PH
   if(q > 0)MAc=Ph0i%*%TH
   Cntc=Ph0i%*%as.matrix(Cnt,k,1)
   zt=da
   #### consider the case t from 1 to pq
   at=matrix((zt[1,]-Cntc),1,k)
   if(pq > 1){
      for (t in 2:pq){
         tmp=matrix((zt[t,]-Cntc),1,k)
         if(p > 0){
            for (j in 1:p){
               if((t-j) > 0){
                  jdx=(j-1)*k
                  tmp1=matrix(zt[(t-j),],1,k)%*%t(as.matrix(ARc[,(jdx+1):(jdx+k)]))
                  tmp=tmp-tmp1
               }
               # end of j-loop
            }
            # end of if(p > 0).
         }
         if(q > 0){
            for (j in 1:q){
               jdx=(j-1)*k
               if((t-j)>0){
                  tmp2=matrix(at[(t-j),],1,k)%*%t(as.matrix(MAc[,(jdx+1):(jdx+k)]))
                  tmp=tmp-tmp2
               }
               #end of j-loop
            }
            #end of if(q > 0)
         }
         at=rbind(at,tmp)
         # end of for(t in 2:pq)
      }
      # end of if(pq > 1) statement
   }
   
   ### for t from ist on
   ist=pq+1
   Pcnt=NULL
   beta=NULL
   if(include.mean){
      beta=matrix(Cntc,1,k)
      Pcnt=c(1)
   }
   if(length(ARc) > 0)beta=rbind(beta,t(ARc))
   if(length(MAc) > 0)beta=rbind(beta,t(MAc))
   idim=k*(p+q)
   if(include.mean){
      Pcnt=c(1)
      idim=idim+1
   }
   #
   for (t in ist:nT){
      Past=NULL
      if(p > 0){
         for (j in 1:p){
            Past=c(Past,zt[(t-j),])
         }
      }
      if( q > 0){
         for (j in 1:q){
            Past=c(Past,at[(t-j),])
         }
      }
      tmp = matrix(c(Pcnt,Past),1,idim)%*%beta
      tmp3=zt[t,]-tmp
      at=rbind(at,tmp3)
   }
   #### skip the first max(p,q) residuals.
   at=as.matrix(at[(ist:nT),])
   sig=crossprod(at,at)/(nT-pq)
   ##
   cat(" ","\n")
   cat("Residuals cov-matrix:","\n")
   print(sig)
   dd=det(sig)
   d1=log(dd)
   ### adjusting for the number of parameters in T-matrix
   jj=0
   for (i in 1:k){
      kdx=c(1:k)[locTmtx[i,] > 1]
      jj=jj+length(kdx)
   }
   aic=d1+2*(npar-jj)/nT
   bic=d1+log(nT)*(npar-jj)/nT
   cat("----","\n")
   cat("aic= ",aic,"\n")
   cat("bic= ",bic,"\n")
   TH1=NULL
   if(length(TH) > 0)TH1=-TH
   
   SCMfit <- list(data=da,SCMs=scms,Tdx=Tdx,locTmtx=locTmtx,locAR=locAR,locMA=locMA,cnst=include.mean,coef=SCMpar,secoef=seSCMpar,residuals=at,Sigma=sig,aic=aic,bic=bic, Ph0=Ph0,Phi=PH,Theta=TH1)
}

"refSCMfit" <- function(model,thres=1.0){
   zt=model$data
   inc.mean=model$cnst
   scms=model$SCMs
   Tdx=model$Tdx
   SCMpar= model$coef
   seSCMpar= model$secoef
   p=max(scms[,1]);q=max(scms[,2])
   #
   nr=length(SCMpar)
   fix=rep(0,nr)
   for (j in 1:nr){
      tt = 0
      if(seSCMpar[j] > 0.000001)tt=SCMpar[j]/seSCMpar[j]
      if(abs(tt) > thres){
         fix[j]=1
      }
      else{
         SCMpar[j]=0
      }
      ### end of "for(j in 1:nr)"
   }
   m1=SCMfit(zt,scms,Tdx,include.mean=inc.mean,fixed=fix,ref=1,SCMpar=SCMpar,seSCMpar=seSCMpar)
   locAR=m1$locAR; locMA=m1$locMA; Tdx=m1$Tdx; locTmtx=m1$locTmtx
   SCMpar=m1$coef; seSCMpar=m1$secoef; scms=m1$SCMs
   sig=m1$Sigma; aic=m1$aic; bic=m1$bic
   Ph0=m1$Ph0; PH=m1$Phi; TH=m1$Theta; if(length(TH)>0)TH=-TH
   at=m1$residuals
   
   refSCMfit <- list(data=zt,SCMs=scms,Tdx=Tdx,locTmtx=locTmtx,locAR=locAR,locMA=locMA,cnst=inc.mean,coef=SCMpar,secoef=seSCMpar,residuals=at,Sigma=sig,aic=aic,bic=bic, Ph0=Ph0,Phi=PH,Theta=TH)
}


################### Co-integration part
"ECMvar1" <- function(x,p,wt,include.const=FALSE,fixed=NULL,output=TRUE){
   # Fits an error-correction VAR model.
   ### This program assumes the co-integrating process w(t) is known.
   ###
   if(!is.matrix(x))x=as.matrix(x)
   if(p < 1)p=1
   nT=dim(x)[1]
   k=dim(x)[2]
   dx=x[2:nT,]-x[1:(nT-1),]
   dx=rbind(rep(0,k),dx)
   wt=as.matrix(wt)
   m=dim(wt)[2]
   wtadj=wt
   ### number of parameters in each equation
   idm=k*(p-1)+m
   if(include.const){
      idm=idm+1
   }
   else{
      wtadj=wt-matrix(1,nT,1)%*%matrix(apply(wt,2,mean),1,m)
   }
   # effective sample size
   ist=max(1,p)
   ne=nT-ist+1
   y=dx[ist:nT,]
   xmtx=wtadj[(ist-1):(nT-1),]
   if(include.const)xmtx=cbind(xmtx,rep(1,(nT-ist+1)))
   if(p > 1){
      for (i in 2:p){
         ii=i-1
         xmtx=cbind(xmtx,dx[(ist-ii):(nT-ii),])
      }
   }
   y=as.matrix(y)
   xmtx=as.matrix(xmtx)
   sdbeta=matrix(1,idm,k)
   #### beta denotes the paramaters in the multiple linear regression format
   #### and is of dimension idm * k
   if(length(fixed) < 1){
      xpx = t(xmtx)%*%xmtx
      xpxinv=solve(xpx)
      xpy=t(xmtx)%*%y
      beta=xpxinv%*%xpy
      yhat=xmtx%*%beta
      resi=y-yhat
      sse=(t(resi)%*%resi)/ne
      dd=diag(xpxinv)
      sdbeta=matrix(1,idm,k)
      for (i in 1:k){
         sdbeta[,i]=sqrt(sse[i,i]*dd)
      }
      npar=idm*k
      ##
   }
   else{
      ## perform estimation equation-by-equation
      resi=NULL
      beta=matrix(0,idm,k)
      sdbeta=matrix(1,idm,k)
      npar=0
      for (i in 1:k){
         idx=c(1:idm)[fixed[,i]==1]
         npi=length(idx)
         cat("Equation: ",i," npar = ",npi,"\n")
         npar=npar+npi
         if(npi > 0){
            xm=xmtx[,idx]
            xpx=t(xm)%*%xm
            xpxinv=solve(xpx)
            dd=diag(xpxinv)
            xpy=t(xm)%*%y[,i]
            betai=xpxinv%*%xpy
            beta[idx,i]=betai
            res=y[,i]-xm%*%betai
            resi=cbind(resi,res)
            se2=sum(res^2)/ne
            sdbeta[idx,i]=sqrt(dd*se2)
            ## end of if(npi > 0)
         }
         # end of for(i in 1:k)
      }
      ##
   }
   sse=(t(resi)%*%resi)/ne
   ### print parameter estimates
   aic=0; bic=0
   if(output){
      alpha=beta[1:m,]
      icnt=m
      if(include.const){
         icnt=m+1
         c=beta[icnt,]
      }
      se=sdbeta[1:m,]
      cat("alpha: ","\n")
      print(t(alpha),digits=3)
      cat("standard error","\n")
      print(t(se),digits=3)
      if(include.const){
         cat("constant term:","\n")
         print(c,digits=3)
         se=sdbeta[icnt,]
         cat("standard error","\n")
         print(se,digits=3)
      }
      ## AR coefficients if any
      if(p > 1){
         cat("AR coefficient matrix","\n")
         jst=icnt
         for (i in 1:(p-1)){
            cat("AR(",i,")-matrix","\n")
            phi=t(beta[(jst+1):(jst+k),])
            se=t(sdbeta[(jst+1):(jst+k),])
            print(phi,digits=3)
            cat("standard error","\n")
            print(se,digits=3)
            jst=jst+k
            ###cat("      ","\n")
         }
         # end of printing AR coefficients
      }
      cat("-----","\n")
      cat("Residuals cov-mtx:","\n")
      print(sse)
      #sse=sse*ne/nT
      cat("      ","\n")
      dd=det(sse)
      cat("det(sse) = ",dd,"\n")
      d1=log(dd)
      aic=d1+(2*npar)/nT
      bic=d1+log(nT)*npar/nT
      cat("AIC = ",aic,"\n")
      cat("BIC = ",bic,"\n")
      ## end of if(output)
   }
   
   ECMvar1 <-list(data=x,wt=wt,arorder=p,include.const=include.const,coef=beta,aic=aic,bic=bic,residuals=resi,secoef=sdbeta,Sigma=sse)
}

###
"refECMvar1" <- function(m1,thres=1.0){
   ### m1 is a fitted model from ECMvar1 or refECMvar1.
   x=m1$data; wt=m1$wt; p=m1$arorder; include.con=m1$include.const
   coef=m1$coef
   secoef=m1$secoef
   idm=dim(coef)[1]
   k=dim(coef)[2]
   fix=matrix(0,idm,k)
   #
   for (i in 1:k){
      tra=coef[,i]/secoef[,i]
      idx=c(1:idm)[abs(tra) > thres]
      fix[idx,i]=1
   }
   #
   mm=ECMvar1(x,p,wt,include.const=include.con,fixed=fix)
   beta=mm$coef; sdbeta=mm$secoef
   aic=mm$aic; bic=mm$bic; resi=mm$residuals; sse=mm$Sigma
   
   refECMvar1 <- list(data=x,wt=wt,arorder=p,include.const=include.con,coef=beta,aic=aic,bic=bic,residuals=resi,secoef=sdbeta,Sigma=sse)
}

#####
"ECMvar" <- function(x,p,ibeta,include.const=FALSE,fixed=NULL,alpha=NULL,se.alpha=NULL,se.beta=NULL,phip=NULL,se.phip=NULL){
   # Fits an error-correction VAR model.
   ### This program assumes the co-integrating process w(t) is unknown.
   ### It is a refined version of ECMvar1.
   ### ibeta: initial estimates of beta-matrix. (k by m matrix). 
   ###   Typically, it is available from the co-integration test.
   if(!is.matrix(x))x=as.matrix(x)
   if(!is.matrix(ibeta))ibeta=as.matrix(ibeta)
   if(p < 1)p=1; m=dim(ibeta)[2]
   cat("Order p: ",p," Co-integrating rank: ",m,"\n")
   nT <- dim(x)[1]; k <- dim(x)[2]
   dx=x[2:nT,]-x[1:(nT-1),]
   dx=rbind(rep(0,k),dx)
   if(length(fixed) < 1){
    #### Obtain initial parameter estimates via ECMvar1, if necessary
    wt=x%*%ibeta
    m1=ECMvar1(x,p,wt,include.const=include.const,output=FALSE)
    est=m1$coef
    se.est=m1$secoef
    alpha=t(est[1:m,])
    se.alpha=t(se.est[1:m,])
    icnt=m
    idm=dim(est)[1]
    ## phip inlcudes the constant term, if any.
    if(idm > icnt){
      phip=est[(icnt+1):idm,]
      se.phip=se.est[(icnt+1):idm,]
      }
    par=c(alpha)
    separ=c(se.alpha)
    par=c(par,c(ibeta[(m+1):k,]))
    separ=c(separ,rep(1/sqrt(nT),(k-m)*m))
    par=c(par,c(phip))
    separ=c(separ,c(se.phip))
    }
   else{
    par=c(alpha); separ=c(se.alpha)
    par=c(par,c(ibeta[(m+1):k,])); separ=c(separ,rep(1/sqrt(nT),(k-m)*m))
    idm=dim(phip)[1]
    for (j in 1:k){
     idx=c(1:idm)[fixed[,j]==1]
     if(length(idx) > 0){
       par=c(par,phip[idx,j]); separ=c(separ,se.phip[idx,j])
       }
      }
    }
   npar=length(par)
##Setup the X-matrix for ECMvar estimation. 
 ECMxmtx <- function(x,p,m,include.const){
   nT <- dim(x)[1]; k <- dim(x)[2]
   dx=x[2:nT,]-x[1:(nT-1),]
   dx=rbind(rep(0,k),dx)
   ist=p
   if(ist < 2)ist=2
   xm=x[(ist-1):(nT-1),]
   ne=nT-ist+1
   if(include.const)xm=cbind(xm,rep(1,ne))
   if(p > 1){
      for (ii in 1:(p-1)){
         xm=cbind(xm,dx[(ist-ii):(nT-ii),])
      }
   }   
   ECMxmtx <- list(xm = xm, y=dx[ist:nT,])
  }
   m2=ECMxmtx(x,p,m,include.const)
   ECMy <- m2$y; ECMxm <- m2$xm
   ##
   cat("Number of parameters: ",length(par),"\n")
   cat("initial estimates: ",par,"\n")
   ### Set up lower and upper bounds
   lowerBounds=par; upperBounds=par
   mult=1.5
   for (j in 1:npar){
      lowerBounds[j] = par[j]-mult*separ[j]
      upperBounds[j] = par[j]+mult*separ[j]
   }
   cat("Par. Lower-bounds: ",lowerBounds,"\n")
   cat("Par. Upper-bounds: ",upperBounds,"\n")

 LECMvar <- function(par,x=x,p=p,m=m,include.const=include.const,fixed=fixed,ECMy=ECMy,ECMxm=ECMxm){
   nT <- dim(x)[1]; k <- dim(x)[2]
   dx=x[2:nT,]-x[1:(nT-1),]
   dx=rbind(rep(0,k),dx)
   #
   km = k*m; kmm = k-m
   npar=length(par)
   if(length(fixed) < 1){
      alpha <- matrix(par[1:km],k,m)
      Im = diag(rep(1,m))
      icnt=(k-m)*m
      beta1=  matrix(par[(km+1):(km+icnt)],kmm,m)
      beta=rbind(Im,beta1)
      Pi=alpha%*%t(beta)
      icnt=icnt+km
      idm=(p-1)*k
      if(include.const)idm=idm+1
      if(npar > icnt){
         ome=matrix(par[(icnt+1):npar],idm,k)
       }
      else{
         ome=NULL
       }
      Ome=rbind(t(Pi),ome)
      resi=ECMy-as.matrix(ECMxm)%*%as.matrix(Ome)
      # end of if(length(fix) < 1)
    }
   else{
      alpha <- matrix(par[1:km],k,m)
      Im <- diag(rep(1,m))
      icnt=(k-m)*m
      beta1 <- matrix(par[(km+1):(km+icnt)],kmm,m)
      beta <- rbind(Im,beta1)
      Pi = alpha%*%t(beta)
      icnt=icnt+km
      idm=(p-1)*k
      if(include.const)idm=idm+1
      if(npar > icnt){
       ome=matrix(0,idm,k)
       for (j in 1:k){
         idx=c(1:idm)[fixed[,j]==1]
         if(length(idx) > 0){
           jj=length(idx)
           ome[idx,j]=par[(icnt+1):(icnt+jj)]
            icnt=icnt+jj
          }
         }
        }
       else{
        ome=NULL
        }
      Ome=rbind(t(Pi),ome)
      resi=ECMy - as.matrix(ECMxm)%*%as.matrix(Ome)
    }
   ## evaluate the log likelihood
   sig=t(resi)%*%resi/nT
   ll=dmvnorm(resi,rep(0,k),sig)
   LECMvar=-sum(log(ll))
   LECMvar
  }
   ###mm=optim(par,LECMvar,method=c("L-BFGS-B"),lower=lowerBounds,upper=upperBounds,hessian=TRUE)
   ###mm=optim(par,LECMvar,method=c("BFGS"),hessian=TRUE)
   ##est=mm$par
   ##H=mm$hessian
   details=FALSE
   # Step 5: Estimate Parameters and Compute Numerically Hessian:
   if(details){ 
      fit = nlminb(start = par, objective = LECMvar,x=x,p=p,m=m,include.const=include.const,fixed=fixed,ECMy=ECMy,ECMxm=ECMxm,
      lower = lowerBounds, upper = upperBounds, control = list(trace=3))
   }
   else {
      fit = nlminb(start = par, objective = LECMvar,x=x,p=p,m=m,include.const=include.const,fixed=fixed,ECMy=ECMy,ECMxm=ECMxm,
       lower = lowerBounds, upper = upperBounds)
   }
   epsilon = 0.0001 * fit$par
   ##npar=length(par)
   Hessian = matrix(0, ncol = npar, nrow = npar)
   for (i in 1:npar) {
      for (j in 1:npar) {
         x1 = x2 = x3 = x4  = fit$par
         x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
         x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
         x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
         x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
         Hessian[i, j] = (LECMvar(x1,x=x,p=p,m=m,include.const=include.const,fixed=fixed,ECMy=ECMy,ECMxm=ECMxm)
                         -LECMvar(x2,x=x,p=p,m=m,include.const=include.const,fixed=fixed,ECMy=ECMy,ECMxm=ECMxm)
                         -LECMvar(x3,x=x,p=p,m=m,include.const=include.const,fixed=fixed,ECMy=ECMy,ECMxm=ECMxm)
                         +LECMvar(x4,x=x,p=p,m=m,include.const=include.const,fixed=fixed,ECMy=ECMy,ECMxm=ECMxm))/
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
   ### print parameter estimates in the ECM model
   km=k*m
   kmm = (k-m)*m
   alpha=matrix(est[1:km],k,m)
   se.alpha=matrix(se.coef[1:km],k,m)
   beta1=matrix(est[(km+1):(km+kmm)],(k-m),m)
   se.beta1=matrix(se.coef[(km+1):(km+kmm)],(k-m),m)
   Im=diag(rep(1,m))
   beta=rbind(Im,beta1); se.beta=rbind(Im,se.beta1)
   icnt=km+kmm
   idm=k*(p-1)
   if(include.const)idm=idm+1
   ### Again, phip includes the constant
   if(icnt < npar){
     if(length(fixed) < 1){
       phip=matrix(est[(icnt+1):npar],idm,k)
       se.phip=matrix(se.coef[(icnt+1):npar],idm,k)
       }
      else{
       phip=matrix(0,idm,k); se.phip=matrix(1,idm,k)
       for (j in 1:k){
        idx=c(1:idm)[fixed[,j]==1]
        jj = length(idx)
        if(jj > 0){
          phip[idx,j]=est[(icnt+1):(icnt+jj)]
          se.phip[idx,j]=se.coef[(icnt+1):(icnt+jj)]
          icnt=icnt+jj
         }
        }
      }
   }   
   cat("alpha: ","\n")
   print(alpha,digits=3)
   cat("standard error","\n")
   print(se.alpha,digits=3)
   cat("beta: ","\n")
   print(beta,digits=3)
   cat("standard error","\n")
   print(se.beta,digits=3)
   #
   icnt=0
   if(include.const){
      cat("constant term:","\n")
      print(phip[1,],digits=3)
      se=se.phip[1,]
      cat("standard error","\n")
      print(se,digits=3)
      icnt=1
   }
   ## AR coefficients if any
   if(p > 1){
      cat("AR coefficient matrix","\n")
      jst=icnt
      for (i in 1:(p-1)){
         cat("AR(",i,")-matrix","\n")
         phi=t(phip[(jst+1):(jst+k),])
         se=t(se.phip[(jst+1):(jst+k),])
         print(phi,digits=3)
         cat("standard error","\n")
         print(se,digits=3)
         jst=jst+k
         ###cat("      ","\n")
      }
      # end of printing AR coefficients
   }
   ## compute the residual covariance matrix
   Pi=alpha%*%t(beta)
   Ome=rbind(t(Pi),phip)
   resi=ECMy-as.matrix(ECMxm)%*%Ome
   sse=t(resi)%*%resi/nT
   #
   cat("-----","\n")
   cat("Residuals cov-mtx:","\n")
   print(sse)
   #sse=sse*ne/nT
   cat("      ","\n")
   dd=det(sse)
   cat("det(sse) = ",dd,"\n")
   d1=log(dd)
   aic=d1+(2*npar)/nT
   bic=d1+log(nT)*npar/nT
   cat("AIC = ",aic,"\n")
   cat("BIC = ",bic,"\n")
   
   ECMvar <-list(data=x,ncoint=m,arorder=p,include.const=include.const,alpha=alpha,se.alpha=se.alpha,beta=beta,se.beta=se.beta,aic=aic,bic=bic,residuals=resi,Phip=phip,se.Phip=se.phip,Sigma=sse)
}



#### refinemeant of ECMvar
"refECMvar" <- function(m1,thres=1.0){
   ### m1 is a fitted model from ECMvar or refECMvar. 
   x=m1$data; m=m1$ncoint; p=m1$arorder; include.const=m1$include.const
   alpha=m1$alpha; se.alpha=m1$se.alpha; beta=m1$beta; se.beta=m1$se.beta
   Phip=m1$Phip; se.Phip=m1$se.Phip; aic=m1$aic; bic=m1$bic; resi=m1$residuals; sse=m1$Sigma
   if(p < 2){
    cat("refinement only applies to the case of p > 1","\n")
    }
   else {
    k <- ncol(x)
    idm=(p-1)*k
    if(include.const){idm=idm+1}
    fixed=matrix(0,idm,k)
    for (j in 1:k){
     tra=Phip[,j]/se.Phip[,j]
     idx=c(1:idm)[abs(tra) > thres]
     fixed[idx,j]=1
     }
    mm=ECMvar(x,p,beta,include.const=include.const,fixed=fixed,alpha=alpha,se.alpha=se.alpha,se.beta=se.beta,phip=Phip,se.phip=se.Phip)
    alpha=mm$alpha; se.alpha=mm$se.alpha; beta=mm$beta; se.beta=mm$se.beta
    Phip=mm$Phip; se.Phip=mm$se.Phip; aic=mm$aic; bic=mm$bic
    resi=mm$residuals; sse=mm$Sigma; include.const=mm$include.const
    }
   refECMvar <-list(data=x,ncoint=m,arorder=p,include.const=include.const,alpha=alpha,se.alpha=se.alpha,beta=beta,se.beta=se.beta,aic=aic,bic=bic,residuals=resi,Phip=Phip,se.Phip=se.Phip,Sigma=sse)
   
}

"SWfore" <- function(y,x,orig,m){
   ### Performs Stock and Watson's diffusion index prediction
   ### y: dependent variable
   ### x: observed regressors
   ### orig: forecast origin
   ### m: selected number of PCs
   ###
   ### Output: Forecasts and MSE of forecasts (if data available)
   if(!is.matrix(x))x=as.matrix(x)
   nT=dim(x)[1]
   k=dim(x)[2]
   if(orig > nT)orig=nT
   if(m > k)m=k; if(m < 1)m=1
   # standardize the predictors
   x1=x[1:orig,]
   me=apply(x1,2,mean)
   se=sqrt(apply(x1,2,var))
   x1=x
   for (i in 1:k){
      x1[,i]=(x1[,i]-me[i])/se[i]
   }
   #
   V1=cov(x1[1:orig,])
   m1=eigen(V1)
   sdev=m1$values
   M=m1$vectors
   M1=M[,1:m]
   Dindex=x1%*%M1
   y1=y[1:orig]; DF=Dindex[1:orig,]
   mm=lm(y1~DF)
   coef=matrix(mm$coefficients,(m+1),1)
   #cat("coefficients: ","\n")
   #print(round(coef,4))
   yhat=NULL; MSE=NULL
   if(orig < nT){
      newx=cbind(rep(1,(nT-orig)),Dindex[(orig+1):nT,])
      yhat=newx%*%coef
      err=y[(orig+1):nT]-yhat
      MSE=mean(err^2)
      cat("MSE of out-of-sample forecasts: ",MSE,"\n")
   }
   
   SWfore <- list(coef=coef,yhat=yhat,MSE=MSE,loadings=M1,DFindex=Dindex)
}

####
"apca" <- function(da,m){
   ### Perform asymptotic PCA when the number of observations is smaller than
   ### the number of variables.
   if(!is.matrix(da))da=as.matrix(da)
   if(m < 1)m=1
   nT=dim(da)[1]
   k=dim(da)[2]
   ### check the validity for performing asymptotic pca.
   if(k <= nT){
      da=t(da)
      nT=dim(da)[1]
      k=dim(da)[2]
   }
   m1=princomp(t(da),cor=F,rotation="none")
   print(summary(m1))
   factors=matrix(m1$loadings,nT,nT)
   factors=factors[,1:m]
   loadings=m1$scores[,1:m]
   sdev=m1$sdev
   
   apca <- list(sdev=sdev,factors=factors,loadings=loadings)
}


#### Constrained factor models of Tsai and Tsay (2011)
####
"hfactor" <- function(X,H,r){
   # Performs estimation of a constrained factor model. The data matrix is "X".
   # The column constraint matrix is H.
   # r: The number of common factor.
   # The program uses a two-step procedure to implement weighted LS estimates.
   # The standardized X follows the following factor model:
   #  SX_t = H * Omega * F_t + epsilon_t
   #
   # This program was written in 2009.
   if(!is.matrix(X))X=as.matrix(X)
   if(!is.matrix(H))H=as.matrix(H)
   N=ncol(X)
   nT=nrow(X)
   m=ncol(H)
   x=X
   print("Data are individually standardized")
   if(r < 1)r=1
   # standardized the data
   x=scale(X,center=TRUE,scale=TRUE)
   V1=cov(x)
   mpca=eigen(V1)
   cat("First r eigenvalues of the correlation matrix: ","\n")
   print(mpca$values[1:r])
   ratio=sum(mpca$values[1:r])/N
   cat("Variability explained: ","\n")
   print(ratio)
   cat("Loadings: ","\n")
   print(mpca$vectors[,1:r],digits=3)
   
   ## New version use square-root of (H'H)^{-1} and 
   ## [(H'H)^{-1/2}H'X'][XH(H'H)^{-1/2}] to perform eigenvalues analysis
   Y=as.matrix(x%*%H)
   HPH=t(H)%*%H
   mh=msqrt(HPH)
   Mhinv=mh$invsqrt
   YH=Y%*%Mhinv
   D=t(YH)%*%YH/nT
   m1=eigen(D)
   cat("eigenvalues of constrained part: ","\n")
   print(m1$values,digits=3)
   d=m1$vectors[,1:r]
   Fhat=YH%*%d
   HPHi=Mhinv%*%Mhinv
   
   # standardize the eigen vectors so that cov(f_t) = I_r.
   V2=apply(Fhat,2,var)
   s1=sqrt(V2)
   for (i in 1:r){
      Fhat[,i]=Fhat[,i]*(1/s1[i])
   }
   
   Omehat=HPHi%*%t(Y)%*%Fhat/nT
   print("Omega-Hat")
   print(Omehat,digits=3)
   
   HO=H%*%Omehat
   pro=sum(diag(t(HO)%*%HO))/N
   cat("Variation explained by the constrained factors: ","\n")
   print(pro)
   cat("H*Omega: constrained loadings ","\n")
   HOOH=HO%*%t(HO)
   Cload=HO
   dd=diag(t(HO)%*%HO)
   for (i in 1:r){
      Cload[,i]=HO[,i]/sqrt(dd[i])
   }
   print(Cload,digits=3)
   Psi=V1-HOOH
   print("Psi:")
   print(Psi,digits=3)
   mpsi=eigen(Psi)
   cat("Diagonal elements of Psi:","\n")
   print(diag(Psi),digits=3)
   cat("eigenvalues of Psi:","\n")
   print(mpsi$values,digits=3)
   
   list(Omega=Omehat,F=Fhat,Psi=Psi)
}

"msqrt" <- function(M){
   # computes the square-root of a positive definite matrix
   if(!is.matrix(M))M=as.matrix(M)
   n1=nrow(M)
   if(n1 == 1){
      Mh=sqrt(M)
      Mhinv=1/Mh
   }
   if(n1 > 1){
      M=(M+t(M))/2
      m1=eigen(M)
      V=m1$vectors
      eiv=sqrt(m1$values)
      L=diag(eiv)
      Linv=diag(1/eiv)
      Mh=V%*%L%*%t(V)
      Mhinv=V%*%Linv%*%t(V)
   }
  msqrt <- list(mtxsqrt=Mh,invsqrt=Mhinv)
}

"BVAR" <- function(z,p=1,C,V0,n0=5,Phi0=NULL,include.mean=T){
## Perform Bayesian estimation of a VAR(p) model
##
## z: time series (T-by-k)
## p: AR order
## phi0: prior mean for coefficient matrix [k-by-(kp+1)]
## C: precision matrix of coefficient matrix. [(kp+1)-by-(kp+1)]
## (V0,n0): prior input for Sigma_a (inverted Wishart parameters)
##
if(!is.matrix(z))z=as.matrix(z)
if(p < 1) p=1
if(!is.matrix(C))C=as.matrix(C)
if(!is.matrix(V0))V0=as.matrix(V0)
if(n0 < 1)n0=1
##
nT=dim(z)[1]
k=dim(z)[2]
idim=k*p+1
if(length(Phi0) <= 0)Phi0=matrix(0,idim,k)
X=NULL
ne=nT-p
if(include.mean)X=rep(1,ne)
for (i in 1:p){
X=cbind(X,z[(p+1-i):(nT-i),])
}
Z=as.matrix(z[(p+1):nT,])
X=as.matrix(X)
### 
XpX=crossprod(X,X)
XpY=crossprod(X,Z)
## Bayesian Estimate
WpW=XpX+C
WpWinv=solve(WpW)
WpY=XpY+C%*%Phi0
Bbhat=WpWinv%*%WpY
bAhat=Z-X%*%Bbhat
bB=Bbhat-Phi0
S=t(bAhat)%*%bAhat +t(bB)%*%C%*%bB
BSig=(V0+S)/(n0+ne-k-1)
SD=kronecker(BSig,WpWinv)
phi=c(Bbhat)
se=sqrt(diag(SD))
Est=cbind(phi,se,phi/se)
colnames(Est) <- c("Est","s.e.","t-ratio")
cat("Bayesian estimate:","\n")
print(Est)
cat("Covariance matrix: ","\n")
print(BSig)
cnst=NULL; Bphi=NULL
if(include.mean){
   cnst=Bbhat[1,]
   Bphi=t(Bbhat[2:idim,])
   }else{
    Bphi=t(Bbhat)
   }
BVAR <- list(phi0=cnst,Phi=Bphi,residuals=bAhat,Sigma=BSig,p=p,priorm=Phi0,precision=C)
}


"comVol" <- function(rtn,m=10,p=1,stand=FALSE){
# checking for common volatility components
if(!is.matrix(rtn))rtn=as.matrix(rtn)
# Fit a VAR(p) model to remove any serial correlations in the data.
if(p < 1){x=scale(rtn,center=T,scale=F)}
 else{
  m1=VAR(rtn,p=p,output=FALSE)
  x=as.matrix(m1$residuals)
 }
#
nT=dim(x)[1]
k=dim(x)[2]
#
if(m < 1)m=1
# standardize the returns
# mean of x is zero because VARfit employs a constant term. 
V1=cov(x)
##print(V1,digits=3)
m1=eigen(V1)
D1=diag(1/sqrt(m1$values))
P1=m1$vectors
Shalf=P1%*%D1%*%t(P1)
x1=x%*%Shalf
#
A=matrix(0,k,k)
for (h in 1:m){
ist=h+1
for (i in 1:k){
for (j in i:k){
Cmtx=matrix(0,k,k)
y2=x1[(ist-h):(nT-h),i]*x1[(ist-h):(nT-h),j]
for (ii in 1:k){
for (jj in ii:k){
y1=x1[ist:nT,ii]*x1[ist:nT,jj]
Cmtx[ii,jj]=cov(y1,y2)*(nT-h)/nT
Cmtx[jj,ii]=Cmtx[ii,jj]
}
}
Cmtx=Cmtx*((nT-h)/nT)
A= A+Cmtx%*%Cmtx
#end of j
}
#end of i
}
#end of h
}
#print(Cmtx)
if(stand){
dd=diag(A)
D=diag(1/sqrt(dd))
A=D%*%A%*%D
}
else{
A=A/(k*(k+1)/2)
}
m2=eigen(A)
Valu=m2$values
Prop=Valu/sum(Valu)
cat("eigen-values: ",Valu,"\n")
cat("proportion:   ",Prop,"\n")
Vec=m2$vectors
Mmtx=Shalf%*%Vec
# normalize each column of Mmtx
for (j in 1:k){
Mmtx[,j]=Mmtx[,j]/sqrt(sum(Mmtx[,j]^2))
}

archTstC <- function(x,m){
# perform F-test for ARCH effect using x^2 series 
# m*Fratio is asymptotically chi-square with m degrees of freedom.
#
if(m < 1)m=1
nT=length(x)
ist=m+1
EffN = nT-m
Xmtx=matrix(1,EffN,1)
Y=matrix(x[ist:nT]^2,EffN,1)
for (j in 1:m){
Xmtx=cbind(Xmtx,x[(ist-j):(nT-j)]^2)
}
XtX=crossprod(Xmtx,Xmtx) 
XtY=crossprod(Xmtx,Y)
beta=solve(XtX,XtY)
Resi=Y-Xmtx%*%beta
Ywm=scale(Y,center=T,scale=F)
SSR=sum(Resi^2)
deg=EffN-m-1
Fratio=((sum(Ywm^2)-SSR)/m)/(SSR/deg)
pv=1-pf(Fratio,m,deg)
result=c(Fratio,pv)
result
}

# Perform ARCH tests for each transformed series
Tst=NULL
Tx = x%*%Mmtx
for (i in 1:k){
TT=NULL
mtst10=archTstC(Tx[,i],10)
TT=c(TT,mtst10)
mtst20=archTstC(Tx[,i],20)
TT=c(TT,mtst20)
mtst30=archTstC(Tx[,i],30)
TT=c(TT,mtst30)
Tst=rbind(Tst,c(i,TT))
}
cat("Checking: ","\n")
cat("Results of individual F-test for ARCH effect","\n")
cat("Numbers of lags used: 10, 20, 30","\n")
cat("Component,(F-ratio P-val) (F-ratio P-val) (F-ratio P-Val)","\n")
print(Tst,digits=3)

comVol <- list(residuals=x,values=m2$values,vectors=m2$vectors,M=Mmtx)
}


"GrangerTest" <- function (X, p = 1, include.mean = T, locInput=c(1)) 
{
    if (!is.matrix(X))X = as.matrix(X)
    Tn = dim(X)[1]
    k = dim(X)[2]
## Re-ordering the components so that the input variables are in front.
    idx=c(1:k)
    if(is.null(locInput))locInput=c(1)
    endog=idx[-locInput]
    jdx=c(locInput,endog)
    x=X[,jdx]
    k1 = length(locInput)
    k2=k-k1
#
    if (p < 1) 
        p = 1
    ne = Tn - p
    ist = p + 1
    y = x[ist:Tn, ]
    if (include.mean) {
        xmtx = cbind(rep(1, ne), x[p:(Tn - 1), ])
    }
    else {
        xmtx = x[p:(Tn - 1), ]
    }
    if (p > 1) {
        for (i in 2:p) {
            xmtx = cbind(xmtx, x[(ist - i):(Tn - i), ])
        }
    }
    ndim = dim(xmtx)[2]
    res = NULL
    xm = as.matrix(xmtx)
    xpx = crossprod(xm, xm)
    xpxinv = solve(xpx)
    xpy = t(xm) %*% as.matrix(y)
    beta = xpxinv %*% xpy
    resi = y - xm %*% beta
    sse = t(resi) %*% resi/(Tn - p - ndim)
    C1 = kronecker(sse, xpxinv)
    bhat = c(beta)
    npar = length(bhat)
    K = NULL
    omega = NULL
##### Locating the zero parameters based on Granger's causality
    for (i in 1:k1) {
      icnt=0
      if(include.mean)icnt=icnt+1
       for (ii in 1:p){
        icnt=icnt+k1
         for (j in 1:k2){
            icnt=icnt+1
            idx=rep(0,npar)
            idx[icnt] = 1
            K = rbind(K, idx)
            omega = c(omega, bhat[icnt])
            }
          }
       }
    K = as.matrix(K)
    v = dim(K)[1]
    cat("Number of targeted zero parameters: ", v, "\n")
    if (v > 0) {
        C2 = K %*% C1 %*% t(K)
        C2inv = solve(C2)
        tmp = C2inv %*% as.matrix(omega, v, 1)
        chi = sum(omega * tmp)
        pvalue = 1 - pchisq(chi, v)
        cat("Chi-square test for Granger Causality and p-value: ", c(chi, pvalue), 
            "\n")
    }
### If p-value is large, perform the estimation of constrained model
 if(pvalue >= 0.05){
  ndim=p*k
  fixed=matrix(1,ndim,k)
  for (i in 1:k1){
    icnt=0
    for (ii in 1:p){
      icnt=icnt+k1
      fixed[(icnt+1):(icnt+k2),i]=0
      }
    }
   if(include.mean){fixed=rbind(rep(0,k),fixed)}
   m1=VAR(X,p=p,include.mean=include.mean,fixed=fixed)
   coef=m1$coef
   secoef=m1$secoef
   aic=m1$aic
   bic=m1$bic
   hq=m1$hq
   resi=m1$residuals
   Sigma=m1$Sigma
   Phi=m1$Phi
   Ph0=m1$Ph0
 }
 else{
   coef=beta
   secoef=NULL
   aic=NULL
   bic=NULL
   hq=NULL
   Sigma=sse
   Phi=NULL
   Ph0=NULL
   }
#
    GrangerTest <- list(data = X, cnst = include.mean, order = p, 
        coef = coef, constraints = K, aic=aic, bic=bic, hq=hq, 
        residuals=resi, secoef=secoef, Sigma=Sigma, 
        Phi=Phi, Ph0=Ph0, omega = omega, covomega = C2, locInput=locInput)
}
