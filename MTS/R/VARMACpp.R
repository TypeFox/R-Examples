"VARMACpp" <- function(da,p=0,q=0,include.mean=T,fixed=NULL,beta=NULL,sebeta=NULL,prelim=F,details=F,thres=2.0){
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
  LLKVARMACpp <- function(par,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed) {
 k    = dim(zt)[2]
nT    = dim(zt)[1]
pqmax = max(p, q)

ListResiduals = .Call("GetVarmaResiduals", zt, fixed, par, p, q, include.mean)		
at  = do.call(rbind, as.list(ListResiduals))
sig = t(at) %*% at/(nT-pqmax)
ll  = dmvnorm(at,rep(0,k),sig)
LLKVARMACpp =-sum(log(ll))   
LLKVARMACpp
}





   # Step 5: Estimate Parameters and Compute Numerically Hessian:
   if(details){
      fit = nlminb(start = par, objective = LLKVARMACpp,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed,
      lower = lowerBounds, upper = upperBounds, control = list(trace=3))
   }
   else {
      fit = nlminb(start = par, objective = LLKVARMACpp, zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed, 
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
         Hessian[i, j] = (LLKVARMACpp(x1,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed)
                         -LLKVARMACpp(x2,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed)
                         -LLKVARMACpp(x3,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed)
                         +LLKVARMACpp(x4,zt=da,p=p,q=q,include.mean=include.mean,fixed=fixed))/
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
   
   VARMACpp <- list(data=da,ARorder=p,MAorder=q,cnst=include.mean,coef=beta,secoef=sebeta,residuals=at,Sigma=sig,aic=aic,bic=bic,Phi=PH,Theta=TH,Ph0=Ph0)
}