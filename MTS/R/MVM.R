#### Multivariate Volatility Modeling
"dccFit" <- function(rt,type="TseTsui",theta=c(0.90,0.02),ub=c(0.92,0.079999),lb=c(0.4,0.0001),cond.dist="std",df=7,m=0){
   # Estimation of DCC model using sample correlation matrix as R0
   # rt: standardized return series
   # theta=c(theta1,theta2): parameter vector
   # ub and lb: upper and lower bounds for theta.
   # R0: initial correlation matrix
   # m = number of past observations used to update the parameters Tse and Tsui model.
   #
   if(!is.matrix(rt))rt=as.matrix(rt)
   nT=dim(rt)[1]; k=dim(rt)[2]
   if (m < 2) m=k+1
   if(length(lb) > 2)lb=lb[1:2]
   if(length(ub) > 2)ub=ub[1:2]
   R0=cor(rt)
   th0=Vech(R0); one1=matrix(1,m,1)
   # Prepare quantities for recursive computing of correlation matrices
   if(type=="TseTsui"){
      loccor=one1%*%th0
      for (t in (m+1):nT){
         v1=cor(rt[(t-m):(t-1),])
         loccor=rbind(loccor,Vech(v1))
      }
   }
   else{
      ### here loccor is just a storage.
      loccor=matrix(th0,1,k*(k+1)/2)
      c1=NULL
      for (i in 1:k){
         c1=cbind(c1,rt[,i:k]*rt[,i])
      }
      loccor=rbind(loccor,c1[-nT,])
   }
   ##
   if(cond.dist=="norm"){
      par=theta
      c1=lb
      c2=ub
   }
   else{
      ### Student-t distribution
      par=c(theta,df)
      c1=c(lb,5.01)
      c2=c(ub,20)
   }
   
   dcclike <- function(par,rt=rt,m=m,type=type,loccor=loccor,cond.dist=cond.dist){
      # evaluate the log-likelihood function of a time series with multivariate
      # Student-t distribution with v degrees of freedom or multivariate normal.
      # par: parameter vectors
      #
      # Other inputs needed:
      nT <- dim(rt)[1]; k <- dim(rt)[2]
      #
      R0 = cor(rt); th0=Vech(R0)
      #
      theta1=par[1];theta2=par[2]; theta0=1-theta1-theta2; wk0=theta0*th0
      loccor=as.matrix(loccor)
      #
      if(type=="TseTsui"){
         ist=m+1; one1=matrix(1,m,1)
         TH=one1%*%th0
         icnt=0; TH1=NULL
         for (i in 1:(k-1)){
            icnt=icnt+1
            TH1=cbind(TH1,rep(1,(nT-m)))
            for (j in (i+1):k){
               icnt=icnt+1
               bi=wk0[icnt]+theta2*loccor[ist:nT,icnt]
               bb=filter(bi,theta1,"r",init=th0[icnt])
               TH1=cbind(TH1,bb)
            }
         }
         TH1=cbind(TH1,rep(1,(nT-m)))
         TH=rbind(TH,TH1)
      }
      else{
         ist=2; k1=k*(k+1)/2
         TH=matrix(th0,1,k1)
         TH1=NULL
         for (i in 1:k1){
            bi=wk0[i]+theta2*loccor[ist:nT,i]
            bb=filter(bi,theta1,"r",init=th0[i])
            TH1=cbind(TH1,bb)
         }
         TH=rbind(TH,TH1)
         ### re-normalization
         for (t in ist:nT){
            Qt=VechM(TH[t,])
            dd=sqrt(diag(Qt))
            #### only lower triangular matrix is used.
            for (i in 1:k){
               for (j in 1:i){
                  Qt[i,j]=Qt[i,j]/(dd[i]*dd[j])
               }
            }
            TH[t,]=Vech(Qt)
         }
         #end of else
      }
      ### Evaluate the log-likelihood function
      llike=0
      #
      if(cond.dist=="norm"){
         for (t in ist:nT){
            sigma1=VechM(TH[t,])
            l1=dmvnorm(rt[t,],mean=rep(0,k),sigma=sigma1,log=TRUE)
            llike=llike-l1
         }
      }
      else{
         ## Multivariate Student-t
         v=par[3]
         k1=log(gamma((v+k)/2))
         k2=log(gamma(v/2))
         k3=k*log((v-2)*pi)/2
         for (t in ist:nT){
            # Find the estimate of theta from the past m data points
            Rt = VechM(TH[t,])
            Rt=as.matrix(Rt)
            yt=matrix(rt[t,],k,1)
            #### The following uses Cholesky decomposition
            b1=chol(Rt); b1=t(b1)
            b1inv=Lminv(b1)
            d1=prod(diag(b1))^2
            yt=matrix(rt[t,],k,1)
            P1=b1inv%*%yt
            tmp1=t(P1)%*%P1
            llike=llike-k1+k2+0.5*log(d1)+k3+(v+k)/2*log(1+tmp1/(v-2))
            #### end of loop for t.
         }
      }
      llike
   }
   
   m1=optim(par,dcclike,method="L-BFGS-B",rt=rt,m=m,type=type,loccor=loccor,cond.dist=cond.dist,lower=c1,upper=c2,hessian=T)
   #m1=optim(par,dcclike,method="BFGS",hessian=T)
   #m1=optim(par,dcclike,method="Nelder-Mead",hessian=T)
   est=m1$par
   H=m1$hessian
   #print(est)
   #print(H)
   Hi=solve(H)
   se=sqrt(diag(Hi))
   cat("Estimates: ",est,"\n")
   cat("st.errors: ",se,"\n")
   cat("t-values:  ",est/se,"\n")
   ##
   dccRho <- function(par,rt=rt,m=m,type=type,loccor=loccor){
      ##
      nT=dim(rt)[1]; k=dim(rt)[2]
      R0 = cor(rt); th0=Vech(R0)
      #
      theta1=par[1];theta2=par[2]; theta0=1-theta1-theta2; wk0=theta0*th0
      loccor=as.matrix(loccor)
      #
      if(type=="TseTsui"){
         ist=m+1; one1=matrix(1,m,1)
         TH=one1%*%th0
         for (t in ist:nT){
            v1 = wk0+theta1*TH[t-1,]+theta2*loccor[t,]
            TH=rbind(TH,v1)
         }
         ## end of if(type=="TseTsui")
      }
      else{
         ist=2
         TH=matrix(th0,1,k*(k+1)/2)
         for (t in ist:nT){
            v1=wk0+theta1*TH[t-1,]+theta2*loccor[t,]
            TH=rbind(TH,v1)
         }
         ### re-normalization
         for (t in ist:nT){
            Qt=VechM(TH[t,])
            dd=sqrt(diag(Qt))
            # Only lower triangular matrices are used.
            for (i in 1:k){
               for (j in 1:i){
                  Qt[i,j]=Qt[i,j]/(dd[i]*dd[j])
               }
            }
            TH[t,]=Vech(Qt)
         }
         #end of else
      }
      ### compute the volatility matrices
      rho.t = NULL
      for (t in 1:nT){
         Rt=VechM(TH[t,])
         rho.t=rbind(rho.t,c(Rt))
      }
      rho.t
   }
   ### Compute the volatility estimates
   rho.t = dccRho(est,rt=rt,m=m,type=type,loccor=loccor)
   
   dccFit <- list(estimates=est,Hessian=H,rho.t=rho.t)
}

"Vech" <- function(mtx){
   # produces half-stacking vector
   # mtx: symmetric
   if(!is.matrix(mtx))mtx=as.matrix(mtx)
   vec=mtx[,1]
   k=nrow(mtx)
   if(k > 1){
      for (j in 2:k){
         vec=c(vec,mtx[j:k,j])
      }
   }
   vec
}

##
"VechM" <- function(vec){
   # from vec to create a symmetric matrix
   m=length(vec)
   k=(-1+sqrt(1+8*m))/2
   mtx=diag(rep(1,k))
   ist=0
   for (j in 1:k){
      mtx[j:k,j]=vec[(ist+1):(ist+k-j+1)]
      ist=ist+(k-j+1)
   }
   #
   for (i in 1:(k-1)){
      for (j in (i+1):k){
         mtx[i,j]=mtx[j,i]
      }
   }
   mtx
}


"Lminv" <- function(L){
   # find the inverse of a lower triangular matrix.
   if(!is.matrix(L))L=as.matrix(L)
   # transform the matrix so that the diagonal elements are 1.
   di=diag(L)
   k=nrow(L)
   # L = L1 * D, where D is diagonal matrix and the diagonal elements of L1 is 1.
   L1=L
   for (i in 1:k){
      L1[i:k,i]=L1[i:k,i]/di[i]
   }
   #
   mtxinv <- function(x){
      # Computes the inverse of a Lower-triangular matrix (with one on the diagonal)
      if(!is.matrix(x))x=as.matrix(x)
      N=dim(x)[2]
      y= diag(N)*2-x
      if(N > 2){
         for (ii in 3:N){
            jj = ii-2
            for (k in 1:jj){
               j=jj-k+1
               tmp=-x[ii,j]
               if ((ii-1)>j){
                  for (i in (j+1):(ii-1)){
                     tmp=tmp-x[ii,i]*y[i,j]
                  }
               }
               y[ii,j]=tmp
            }
         }
      }
      y
   }
   #
   Linv=mtxinv(L1)
   for (i in 1:k){
      Linv[i,1:i]=Linv[i,1:i]/di[i]
   }
   Linv
}

##########
"dccPre" <- function(rtn,include.mean=T,p=0,cond.dist="norm"){
   ## Fits marginal GARCH models to each component to obtain
   ### marginally standardized series for "dccFit".
   if(!is.matrix(rtn))rtn=as.matrix(rtn)
   RTN=rtn
   if(p > 0){
      m1=VAR(rtn,p)
      RTN=m1$residuals
   }
   else{
      if(include.mean){
         mu=apply(rtn,2,mean)
         RTN=scale(rtn,center=T,scale=F)
         cat("Sample mean of the returns: ",mu,"\n")
         RTN=as.matrix(RTN)
      }
   }
   k=dim(RTN)[2]
   Vol=NULL; sresi=NULL
   Mtxest=NULL; Mtxse=NULL; Mtxtval=NULL
   for (i in 1:k){
      m2=garchFit(~garch(1,1),data=RTN[,i],include.mean=F,trace=F,cond.dist=cond.dist)
      Mtxest=rbind(Mtxest,m2@fit$par); Mtxse=rbind(Mtxse,m2@fit$se.coef)
      Mtxtval=rbind(Mtxtval,m2@fit$tval)
      Vol=cbind(Vol,m2@sigma.t)
      sresi=cbind(sresi,m2@residuals/m2@sigma.t)
      cat("Component: ",i,"\n")
      cat("Estimates: ",round(Mtxest[i,],6),"\n")
      cat("se.coef  : ",round(Mtxse[i,],6),"\n")
      cat("t-value  : ",round(Mtxtval[i,],6),"\n")
   }
   dccPre <- list(marVol=Vol,sresi=sresi,est=Mtxest,se.coef=Mtxse)
}

#############
"EWMAvol" <- function(rtn,lambda=0.96){
   ## Compute exponentially weighted moving average covariance matrix.
   ## If lambda <= 0, likelihood function is used to estimate lambda.
   if(!is.matrix(rtn))rtn=as.matrix(rtn)
   nT=dim(rtn)[1]
   k=dim(rtn)[2]
   x=scale(rtn, center=TRUE,scale=FALSE)
   #
   par=lambda
   MGAUS <- function(par,x=x){
      lambda=par[1]
      h1=1-lambda
      Sigt=cov(x)
      lk=0
      nT=dim(x)[1]
      k=dim(x)[2]
      for (t in 2:nT){
         xx=as.numeric(x[t-1,])
         for (i in 1:k){
            Sigt[i,]=h1*xx[i]*xx+lambda*Sigt[i,]
         }
         ll=dmvnorm(x[t,],rep(0,k),sigma=Sigt,log=TRUE)
         lk=lk-ll
      }
      lk
   }
   if(lambda <= 0){
      # perform QMLE of lambda
      par=c(lambda=0.96)
      S=10^{-5}
      lowerb=c(lambda = S); upperb=c(lambda=1-S)
      fit = nlminb(start = par, objective = MGAUS,x=x,
      lower = lowerb, upper = upperb)
      #
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
            Hessian[i, j] = (MGAUS(x1,x=x)-MGAUS(x2,x=x)-MGAUS(x3,x=x)+MGAUS(x4,x=x))/
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
      lambda=fit$par[1]
   }
   ##
   if(lambda > 0){
      h1=1-lambda
      Sigt=cov(x)
      V1=c(Sigt)
      for (t in 2:nT){
         xx=as.numeric(x[t-1,])
         for (i in 1:k){
            Sigt[i,]= h1*xx*xx[i]+lambda*Sigt[i,]
         }
         V1=rbind(V1,c(Sigt))
      }
   }
   #
   
   EWMAvol <- list(Sigma.t=V1,return=rtn,lambda=lambda)
}

#################
"MCHdiag" <- function(at,Sigma.t,m=10){
   ### Perform diagnostic checks for a fitted volatility models.
   if(!is.matrix(at))at=as.matrix(at)
   if(!is.matrix(Sigma.t))Sigma.t=as.matrix(Sigma.t)
   ##
   nT=dim(at)[1]; k=dim(at)[2]
   nT1=dim(Sigma.t)[1]; k1=dim(Sigma.t)[2]
   if((nT != nT1) || (k1 != k^2)){
      cat("Inconsistency in dimensions","\n")
   }
   else{
      ## perform model checking
      et=NULL
      etat=NULL
      for (i in 1:nT){
         Vt=matrix(Sigma.t[i,],k,k)
         Vtinv=solve(Vt)
         x=matrix(at[i,],1,k)
         tmp=x%*%Vtinv%*%t(x)-k
         et=c(et,tmp)
         m1=eigen(Vt)
         P=m1$vectors; lam=m1$values
         d1=diag(1/sqrt(lam))
         Vthalf=P%*%d1%*%t(P)
         wk=x%*%Vthalf
         etat=rbind(etat,wk)
      }
      m1=acf(et,m,plot=FALSE)
      acf=m1$acf[2:(m+1)]
      tmp=acf^2/c(rep(nT,m)-c(1:m))
      Q1=sum(tmp)*nT*(nT+2)
      pv1=1-pchisq(Q1,m)
      #### Prepare for rank-based test
      lag=m
      mu=-(rep(nT,lag)-c(1:lag))/(nT*(nT-1))
      v1=rep(5*nT^4,lag)-(5*c(1:lag)+9)*nT^3+9*(c(1:lag)-2)*nT^2+2*c(1:lag)*(5*c(1:lag)+8)*nT+16*c(1:lag)^2
      v1=v1/(5*(nT-1)^2*nT^2*(nT+1))
      ret=rank(et)
      m2=acf(ret,m,plot=FALSE)
      acf=m2$acf[2:(m+1)]
      Qr=sum((acf-mu)^2/v1)
      pv2=1-pchisq(Qr,m)
      ###
      ### Multivariate portmanteau test
      x=etat^2
      g0=var(x)
      ginv=solve(g0)
      qm=0.0
      for (i in 1:lag){
         x1=x[(i+1):nT,]
         x2=x[1:(nT-i),]
         g = cov(x1,x2)
         g = g*(nT-i-1)/(nT-1)
         h=t(g)%*%ginv%*%g%*%ginv
         qm=qm+nT*nT*sum(diag(h))/(nT-i)
      }
      QKm=qm; pv3=1-pchisq(QKm,k^2*m)
      ### Robust multivariate tests
      q95=quantile(et,0.95)
      idx=c(1:nT)[et <= q95]
      x=etat[idx,]^2
      eT=length(idx)
      g0=var(x)
      ginv=solve(g0)
      qm=0.0
      for (i in 1:lag){
         x1=x[(i+1):eT,]
         x2=x[1:(eT-i),]
         g = cov(x1,x2)
         g = g*(eT-i-1)/(eT-1)
         h=t(g)%*%ginv%*%g%*%ginv
         qm=qm+eT*eT*sum(diag(h))/(eT-i)
      }
      Qrm = qm; pv4=1-pchisq(Qrm,k^2*m)
      ##
      cat("Test results: ","\n")
      cat("Q(m) of et:","\n")
      cat("Test and p-value: ",c(Q1,pv1),"\n")
      cat("Rank-based test:","\n")
      cat("Test and p-value: ",c(Qr,pv2),"\n")
      cat("Qk(m) of epsilon_t:","\n")
      cat("Test and p-value: ",c(QKm,pv3),"\n")
      cat("Robust Qk(m): ","\n")
      cat("Test and p-value: ",c(Qrm,pv4),"\n")
   }
   ## end of the program
}

#############################
"MCholV" <- function(rtn,size=36,lambda=0.96,p=0){
   ### Perform multivariate volatility modeling via Cholesky Decomposition.
   ### It makes use of the sequential nature of the decomposition.
   ### This is a rather highly structured multivariate volatility model by
   ### keeping the number of parameters low.
   ####
   if(!is.matrix(rtn))rtn=as.matrix(rtn)
   if(size <= 0)size=36
   ### aimed for monthly returns with a window size of 3 years.
   if(p < 0)p=0
   mu=apply(rtn,2,mean)
   if(p==0){
      RTN=scale(rtn,center=T,scale=F)
      cat("Sample means: ",round(mu,6),"\n")
   }
   else{
      m1=VAR(rtn,p)
      RTN=m1$residuals
   }
   nT=dim(RTN)[1]; k=dim(RTN)[2]
   ### Perform volatility estimation: Make use of the sequential nature of the decomposition
   Mtxcoef=matrix(0,k,3); Mtxse=matrix(0,k,3); Mtxtval=matrix(0,k,3)
   bt=RTN[(size+1):nT,1]
   m1=garchFit(~garch(1,1),data=bt,include.mean=F,trace=F)
   Mtxcoef[1,]=m1@fit$par; Mtxse[1,]=c(m1@fit$se.coef); Mtxtval[1,]=m1@fit$tval
   cat("Estimation of the first component","\n")
   cat("Estimate (alpha0, alpha1, beta1): ",round(Mtxcoef[1,],6),"\n")
   cat("s.e.                            : ",round(Mtxse[1,],6),"\n")
   cat("t-value                         : ",round(Mtxtval[1,],6),"\n")
   #
   v1=m1@sigma.t; VOL=v1*v1
   ### BetaU stores the smoothed beta_{ij,t} used.
   BetaU=NULL
   for (i in 2:k){
      ### Perform recursive least squares estimation to obtain beta_{ij,t}.
      betat=NULL
      y=RTN[,i]; x=RTN[,1:(i-1)]
      m2=RLS(y,x,ist=size)
      betat=cbind(betat,-m2$beta)
      bt=cbind(bt,m2$resi)
      ### Perform Volatility estimation of the i-th component
      ##### First, obtain smoothed beta_{ij,t}
      betat=as.matrix(betat)
      BetaM=apply(betat,2,mean)
      betat=scale(betat,center=T,scale=F)
      RTN1 <- RTN[(size+1):nT,1:i]
      k1=dim(betat)[2]; nT1=dim(betat)[1]
      Beta=NULL
      for (j in 1:k1){
         b1=c(betat[1,j],(1-lambda)*betat[-nT1,j])
         bb=filter(b1,lambda,"r",init=0)+BetaM[j]
         Beta=cbind(Beta,bb)
      }
      ### Obtain the residuals
      BetaU=cbind(BetaU,Beta) # store the smoothed beta_{ij,t}
      wk=cbind(Beta,rep(1,nT1))
      bit=apply(RTN1*wk,1,sum)
      bt[,i]=bit  # replace the i-th residuals with smoothed-beta residuals
      m3=garchFit(~garch(1,1),data=bit,include.mean=F,trace=F)
      cat("Component ",i," Estimation Results (residual series):","\n")
      Mtxcoef[i,]=m3@fit$par; Mtxse[i,]=c(m3@fit$se.coef); Mtxtval[i,]=m3@fit$tval
      cat("Estimate (alpha0, alpha1, beta1): ",round(Mtxcoef[i,],6),"\n")
      cat("s.e.                            : ",round(Mtxse[i,],6),"\n")
      cat("t-value                         : ",round(Mtxtval[i,],6),"\n")
      VOL=cbind(VOL,(m3@sigma.t)^2)
   }
   BetaU=as.matrix(BetaU)
   
   LinV <- function(A){
      ### A: lower triangular with one on the diagonal
      ### B: inverse(A)
      if(!is.matrix(A))A=as.matrix(A)
      k1=dim(A)[1]; k2=dim(A)[2]
      k=min(k1,k2)
      for (i in 1:k){
         A[i,i]=1
      }
      B=diag(k)
      for (i in 2:k){
         B[i,i-1]=-A[i,i-1]
      }
      if(k > 2){
         for (i in 3:k){
            iend=i-1  # number of completed band.
            for (j in 1:(k-iend)){
               tmp=0
               for (ii in j:(j+iend-1)){
                  tmp=tmp+A[i+j-1,ii]*B[ii,j]
               }
               B[i+j-1,j]=-tmp
            }# end of j-loop
         } #end of i-loop
      }  # end of (k > 2)
      B
   }
   
   ### Compute the Volatility matrices
   Sigma.t=NULL
   for (it in 1:nT1){
      A=diag(k)
      icnt=0
      for (j in 2:k){
         A[j,1:(j-1)]=BetaU[it,(icnt+1):(icnt+j-1)]
         icnt=icnt+j-1
      }
      B=LinV(A)
      for (i in 1:k){
         A[,i]=VOL[it,i]*B[,i]
      }
      Sigma.t = rbind(Sigma.t,c(A%*%t(B)))
   }
   
   MCholV <- list(betat=-BetaU,bt=bt,Vol=VOL,Sigma.t=Sigma.t)
}


"RLS" <- function(y,x,ist=30,xpxi=NULL,xpy0=NULL){
   ## Compute the recursive least squares estimates of the multiple linear
   ## regression: y=beta'x+e.
   ######## x should include a column of 1's if constant is needed.
   #### xpxi and xpy0 are initial inverse(X'X) and X'y matrices if available.
   #### Return the time-varying least squares estimates and the residual series.
   ####
   if(!is.matrix(x))x=as.matrix(x)
   nT=dim(x)[1]
   k=dim(x)[2]
   #
   if(is.null(xpxi)){
      if(ist <= k)ist=k+1
      xpx0=t(x[1:ist,])%*%x[1:ist,]
      xpxi=solve(xpx0)
      xpy0=t(x[1:ist,])%*%matrix(y[1:ist],ist,1)
   }
   ist=ist+1
   resi=NULL
   beta=matrix(c(xpxi%*%xpy0),1,k)
   for (t in ist:nT){
      jdx=t-ist+1
      x1=matrix(x[t,],1,k)
      wk=x1%*%xpxi
      k1=1+sum(wk*x1)
      wk1=t(wk)%*%wk/k1
      xpxi=xpxi-wk1
      Kmtx=xpxi%*%t(x1)
      tmp=beta[jdx,]-t(Kmtx)*(sum(x1*beta[jdx,])-y[t])
      beta=rbind(beta,tmp)
      res=y[t]-tmp%*%t(x1)
      resi=c(resi,res)
   }
   beta=beta[-1,]
   
   RLS <- list(beta=beta,resi=resi)
}

##################################
"mtCopula" <- function(rt,g1,g2,grp=NULL,th0=NULL,m=0,include.th0=TRUE){
   # Estimation of t-copula when the correlation matrix is constrained.
   # rt: standardized return series
   # grp: vector of group sizes.
   # g1 and g2: initial parameter values
   # th0: initial values of the angles.
   ### Modified in January 2013 for the book. Remove grouping temporarily.
   #
   if(!is.matrix(rt))rt=as.matrix(rt)
   k <- dim(rt)[2]; nT <- dim(rt)[1]
   if(is.null(grp))grp=rep(1,k)
   ### obtain initial value of th0 if needed.
   ### Multiple subroutines are listed next.
   #
   CorTheta <- function(R){
      # Obtain the matrix of angles that give rise to the given correlation matrix
      R1 = R
      if(!is.matrix(R1))R1=as.matrix(R1)
      k=nrow(R1)
      mm=chol(R1)
      R1=as.matrix(mm)
      #
      if(k > 1){
         km1=k-1
         TH=matrix(0,km1,km1)
         for (ii in 1:km1){
            for (j in (ii+1):k){
               TH[ii,j-1]=acos(R1[ii,j])
               for (i in (ii+1):j){
                  R1[i,j]=R1[i,j]/sin(TH[ii,j-1])
               }#end of i
            }#end of j
         }#end of ii
         # vectorization of theta
         theta=TH[1,]
         if(km1 > 1){
            for (i in 2:km1){
               theta=c(theta,TH[i,i:km1])
            }
         }
      }#end of the if statment of (k > 1)
      CorTheta <- list(dim=k, Theta = theta, THmtx=TH)
   }
   
   vtimeU <- function(vec,Up){
      # multiplication of a row vector by an upper triangular matrix
      if(!is.matrix(Up))Up=as.matrix(Up)
      Resu=NULL
      n1=length(vec)
      n2=nrow(Up)
      if(n1==n2){
         Resu=vec[1]*Up[1,1]
         for (i in 2:n1){
            Resu=c(Resu,sum(vec[1:i]*Up[1:i,i]))
         }
         #end of if(n1==n2)
      }
      Resu
   }
   
   TH2idth <- function(TH,grp){
      # Given the TH-matrix and group size, obtains the independent theta.
      #
      if(!is.matrix(TH))TH=as.matrix(TH)
      th=NULL
      k=sum(grp)
      g=length(grp)
      irow = 1
      icnt=0
      for (i in 1:g){
         jcnt=0
         if(i > 1)jcnt=sum(grp[1:(i-1)])
         if(grp[i]>1){
            icnt=icnt+1
            th=c(th,TH[irow,jcnt+1])
            jcnt=jcnt+grp[i]
         }
         if(i < g){
            for (j in (i+1):g){
               icnt=icnt+1
               th=c(th,TH[irow,jcnt+1])
               jcnt=jcnt+grp[j]
            }
         }
         irow=irow+grp[i]
      }#end of i
      th
   }
   
   the2Xmtx <- function(theta,grp){
      # For a given grp and theta, computes the Cholesky's decomposition of the
      # correlation matrix.
      #
      # grp: group sizes.
      #   Data should be in the proper column ordering.
      #
      # Quantites must be computed: Cmtx (correlation matrix),Xmtx (upper triangular)
      # TH (Theta matrix).
      #
      g=length(grp)
      k=sum(grp)
      km1=k-1
      # g: number of groups
      # k: dimension
      # km1: dimension of TH-matrix
      nth=g*(g-1)/2
      for (i in 1:g){
         if(grp[i]>1)nth=nth+1
      }
      # nth: length of the theta-vector
      if(nth != length(theta))print(c(nth,length(theta)))
      Xmtx=matrix(1,k,k)
      for (j in 1:(k-1)){
         for (i in (j+1):k){
            Xmtx[i,j]=0
         }
      }
      Cmtx=t(Xmtx)
      THini=4
      TH=matrix(THini,km1,km1)
      ## TH: upper triangular matrix should be sufficient, but do not
      ## consider the space issue in this program.
      ## Xmtx: In the program, it is a symmetric matrix for ease in computing.
      # icnt: counting the number of elements in theta is used.
      icnt=0
      ## Take the approach: theta ==> TH(independent part) ==> fill TH
      idx = 0
      for (i in 1:g){
         idx=idx+1
         ni=grp[i]
         if(ni > 1){
            icnt=icnt+1
            TH[idx,idx]=theta[icnt]
            # end if (if(ni > 1)
         }
         if(i < g){
            for (j in (i+1):g){
               jdx=sum(grp[1:(j-1)])
               icnt=icnt+1
               TH[idx,jdx]=theta[icnt]
               # end of j
            }
            # end of if(i < g)
         }
         idx=sum(grp[1:i])
         # end of i
      }
      ## Fill in missing theta in the TH matrix
      ### and construct the Xmtx in the process
      #### First, complete the first row
      if(km1 > 1){
         for (j in 2:km1){
            if(TH[1,j]==THini)TH[1,j]=TH[1,j-1]
         }
      }
      for (j in 2:k){
         Xmtx[1,j]=cos(TH[1,j-1])
         Cmtx[j,1]=Xmtx[1,j]
      }
      if(k >1){
         for (j in 2:k){
            for (i in 2:j){
               Xmtx[i,j]=Xmtx[i,j]*sin(TH[1,j-1])
               Cmtx[j,i]=Xmtx[i,j]
            }
         }
      }
      ### complete TH, starting wirh row 2:
      if(km1 > 1){
         for (i in 2:km1){
            for (j in i:km1){
               if(TH[i,j]==THini){
                  if((TH[i-1,j]==TH[i-1,j-1])&&(i<j)){
                     TH[i,j]=TH[i,j-1]
                  }
                  else{
                     # need to fill in: use coef at [i,j]= coef at [i-1,j]
                     c1=0.0
                     for (jj in 1:(i-1)){
                        c1=c1+Cmtx[i-1,jj]*Xmtx[jj,j+1]
                     }
                     for (ii in 1:(i-1)){
                        c1=c1-Cmtx[i,ii]*Xmtx[ii,j+1]
                     }
                     TH[i,j]=acos(c1/(Cmtx[i,i]*Xmtx[i,j+1]))
                  }
               }
               #end of j
            }
            ### construct Xmtx for the i-th row
            for (j in (i+1):k){
               Xmtx[i,j]=Xmtx[i,j]*cos(TH[i,j-1])
               Cmtx[j,i]=Xmtx[i,j]
            }
            for (j in (i+1):k){
               for (ii in (i+1):j){
                  Xmtx[ii,j]=Xmtx[ii,j]*sin(TH[i,j-1])
                  Cmtx[j,ii]=Xmtx[ii,j]
               }
            }
            # end of i
         }
         #end of if(km1 > 1)
      }
      #end of the program
      tmp=1.0
      for (j in 1:(k-1)){
         tmp=tmp*sin(TH[j,km1])
      }
      Xmtx[k,k]=tmp
      Xmtx
   }
   
   
   if(length(th0)==0){
      mi=SCCor(rt,nT,nT,grp)$conCor
      THmi=CorTheta(mi)$THmtx
      th0=TH2idth(THmi,grp)
   }
   if(!include.th0){
      cat("Value of angles: ","\n")
      print(th0)
   }
   #
   g=length(grp)
   ncor=g*(g-1)/2
   for (i in 1:g){
      if(grp[i]>1)ncor=ncor+1
   }
   
   if(m <= 0)m=ncor+1
   ist=m+1
   mgsize=max(grp)
   Theta = matrix(1,m,1)%*%matrix(th0,1,length(th0))
   ### Compute the local information of theta
   for (it in ist:nT){
      # Find the estimate of theta from the past m data points
      V1=SCCor(rt,it-1,m,grp)$conCor
      Tmp=CorTheta(V1)$THmtx
      if(mgsize==1){
         thet=Vech(t(Tmp))
      }
      else{
         thet=TH2idth(Tmp,grp)
      }
      Theta=rbind(Theta,thet)
   }
   #
   if(include.th0){
      par=c(7.0,g1,g2,th0)
      c1=c(5.1,0.2,0.0001,th0*0.8)
      c2=c(20,.95,.04999999,th0*1.1)
   }
   else{
      par=c(7.0,g1,g2)
      c1=c(5.1,0.2,0.00001)
      c2=c(20,.95,0.0499999)
   }
   cat("Lower limits: ",c1,"\n")
   cat("Upper limits: ",c2,"\n")
   
   ### Likelihood function
   
   "mtlikeC" <- function(par,rt=rt,grp=grp,ncor=ncor,m=m,th0=th0,include.th0=include.th0){
      # evaluate the log-likelihood function of a t-copula
      # Student-t distribution with degrees of freedom v.
      # par: parameter vectors, i.e. par=c(v,g1,g2)
      #
      # Other inputs needed:
      # rt: standardized return series nT-by-k
      # th0: initial value of the independent theta-vector
      # grp: grouping of returns
      #
      nT=dim(rt)[1]; k=dim(rt)[2]
      #
      if(include.th0){
         Th0=matrix(par[4:length(par)],1,length(par)-3)
      }
      else{
         Th0 =matrix(th0,1,length(th0))
      }
      on1=matrix(1,m,1)
      TH=on1%*%Th0
      ## ist for the starting point to evaluate the log likelihood function
      ist=m+1
      #
      llike=0
      # Theta: is the estimate of theta from the most recent m data points
      k1=lgamma((par[1]+k)/2)
      k2=lgamma(par[1]/2)
      k4=lgamma((par[1]+1)/2)
      #
      if(include.th0){
         th1=par[4:length(par)]
      }
      else{
         th1 <- th0
      }
      par0=1-par[2]-par[3]; tmp0=par0*th1
      #
      for (t in ist:nT){
         thet=Theta[t,]
         tht=tmp0+TH[t-1,]*par[2]+par[3]*thet
         TH=rbind(TH,tht)
         # Xmtx is upper triangular
         Xmtx=the2Xmtx(tht,grp)
         # d1^2 is the determinant of R_t matrix (|R_t| = |Xmtx'|*|Xmtx|)
         # d1 is the square-root of the determinant
         d1=abs(prod(diag(Xmtx)))
         Xtinv=t(Lminv(t(Xmtx)))
         yt=rt[t,]
         tmp1=vtimeU(yt,Xtinv)
         tmp=sum(tmp1^2)
         #
         llike=llike+k1+(k-1)*k2-k*k4-log(d1)-((par[1]+k)/2)*log(1+tmp/(par[1]-2))
         #
         llike=llike+((par[1]+1)/2)*sum(log(1+yt^2/(par[1]-2)))
         # end of loop for t.
      }
      llike = -llike
      
      llike
   }#end of the program
   
   m1=optim(par,mtlikeC,method="L-BFGS-B",rt=rt,grp=grp,ncor=ncor,m=m,th0=th0,include.th0=include.th0,lower=c1,upper=c2,hessian=T)
   est=m1$par
   H=m1$hessian
   Hi=solve(H)
   se.coef=sqrt(diag(Hi))
   cat("estimates:  ",est,"\n")
   cat("std.errors: ",se.coef,"\n")
   cat("t-values:   ",est/se.coef,"\n")
   ### Alternative approach
   est=m1$par
   npar=length(par)
   epsilon=0.0001*est; epsilon[1]=0.3*est[1]
   Hessian = matrix(0, ncol = npar, nrow = npar)
   for (i in 1:npar) {
      for (j in 1:npar) {
         x1 = x2 = x3 = x4  = m1$par
         x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
         x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
         x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
         x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
         Hessian[i, j] = (mtlikeC(x1,rt=rt,grp=grp,ncor=ncor,m=m,th0=th0,include.th0=include.th0)
         -mtlikeC(x2,rt=rt,grp=grp,ncor=ncor,m=m,th0=th0,include.th0=include.th0)
         -mtlikeC(x3,rt=rt,grp=grp,ncor=ncor,m=m,th0=th0,include.th0=include.th0)
         +mtlikeC(x4,rt=rt,grp=grp,ncor=ncor,m=m,th0=th0,include.th0=include.th0))/
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
   cat("Alternative numerical estimates of se:", "\n")
   cat("st.errors: ",se.coef,"\n")
   cat("t-values:  ",est/se.coef,"\n")
   ### program to compute the time-varying rho.t based on the final estimate
   "mtCopulaVol" <- function(par,include.th0,th0,Theta,ncor,grp,m){
      ##
      if(include.th0){
         th1=par[4:length(par)]
      }
      else{
         th1 <- th0
      }
      par0=1-par[2]-par[3]; tmp0=par0*th1
      ist=m+1
      nT=dim(Theta)[1]
      #
      Xmtx = the2Xmtx(th1,grp)
      Rt=t(Xmtx)%*%Xmtx
      rho.t=NULL
      TH=NULL
      for (i in 1:m){
         rho.t=rbind(rho.t,c(Rt))
         TH=rbind(TH,th1)
      }
      #
      for (t in ist:nT){
         thet=Theta[t,]
         tht=tmp0+TH[t-1,]*par[2]+par[3]*thet
         TH=rbind(TH,tht)
         # Xmtx is upper triangular
         Xmtx=the2Xmtx(tht,grp)
         Rt=t(Xmtx)%*%Xmtx
         rho.t=rbind(rho.t,c(Rt))
      }
      mtCopulaVol <- list(rho.t=rho.t,angles=TH)
   }
   
   mf = mtCopulaVol(est,include.th0,th0,Theta,ncor,grp,m)
   rho.t=mf$rho.t; angles=mf$angles
   
   mtCopula <- list(estimates=est,Hessian=H,rho.t=rho.t,theta.t=angles)
}

###
"SCCor" <- function(rt,end,span,grp){
   # estimation of constrained sample correlation matrix
   # rt: time series T-by-k
   # end: end time point of the window
   # span: window length used to estimate the correlation (span >= k)
   # grp: groups of time series, e.g., grp=c(3,2) indicates the first group consists
   #     of the first three series and the second group the last 2 series.
   if(!is.matrix(rt))rt=as.matrix(rt)
   nT=dim(rt)[1]; k=dim(rt)[2]
   if(end < span)end=nT
   start=end-span+1
   maxgsize=max(grp)
   #
   if(maxgsize > 1){
      #### case of group constraints
      ng=length(grp)
      ncor=ng*(ng-1)/2
      for (i in 1:ng){
         if(grp[i]>1)ncor=ncor+1
      }
      if(span < (ncor+1))span=ncor+1
      if(end < span)end=nT
      start=end-span+1
      V1=cor(rt[start:end,])
      V2=V1
      
      for (i in 1:ng){
         ist=0
         if(i > 1)ist=sum(grp[1:(i-1)])
         #
         if(grp[i]>1){
            ## for j = i case
            tot=0
            icnt=0
            for (ii in (ist+1):(ist+grp[i]-1)){
               for (jj in (ii+1):(ist+grp[i])){
                  icnt=icnt+1
                  tot=tot+V1[ii,jj]
               }
            }
            av=tot/icnt
            for (ii in (ist+1):(ist+grp[i]-1)){
               for (jj in (ii+1):(ist+grp[i])){
                  V2[ii,jj]=av
                  V2[jj,ii]=av
               }
            }
            # end of if grp[i] > 1.
         }
         # for i .ne. j case.
         if((i+1) <= ng){
            for (j in (i+1):ng){
               jst=sum(grp[1:(j-1)])
               tot=0
               icnt=0
               for (ii in (ist+1):(ist+grp[i])){
                  for (jj in (jst+1):(jst+grp[j])){
                     tot=tot+V1[ii,jj]
                     icnt=icnt+1
                  }
               }
               av=tot/icnt
               for (ii in (ist+1):(ist+grp[i])){
                  for (jj in (jst+1):(jst+grp[j])){
                     V2[ii,jj]=av
                     V2[jj,ii]=av
                  }
               }
            }#end of j
         }# end of if statement (i+1) <= ng.
      }
   }
   else{
      ## no constraints
      V1=cor(rt[start:end,])
      V2=V1
   }
   SCCor <- list(unconCor=V1,conCor=V2)
} ### end of SCCor routine



###########################
"archTest" <- function(rt,lag=10){
   ### Perfrom test to detect ARCH effects in a univariate time series rt.
   #### Two types of test are performed.
   #### (a) McLeod-Li (or equivalent Engle's LM) test of rt^2
   #### (b) Rank-based test of Dufour and Roy (1986), uses ranks of rt^2.
   ####
   at=rt
   if(is.matrix(at))at=at[,1]
   m1=acf(at^2,lag.max=lag,plot=F)
   acf=m1$acf[2:(lag+1)]
   nT=length(at)
   c1=c(1:lag)
   deno=rep(nT,lag)-c1
   Q=sum(acf^2/deno)*nT*(nT+2)
   pv1=1-pchisq(Q,lag)
   ###
   cat("Q(m) of squared series(LM test): ","\n")
   cat("Test statistic: ",Q," p-value: ",pv1,"\n")
   ##
   rk=rank(at^2)
   m2=acf(rk,lag.max=lag,plot=F)
   acf=m2$acf[2:(lag+1)]
   mu=-(rep(nT,lag)-c(1:lag))/(nT*(nT-1))
   v1=rep(5*nT^4,lag)-(5*c(1:lag)+9)*nT^3+9*(c(1:lag)-2)*nT^2+2*c(1:lag)*(5*c(1:lag)+8)*nT+16*c(1:lag)^2
   v1=v1/(5*(nT-1)^2*nT^2*(nT+1))
   QR=sum((acf-mu)^2/v1)
   pv2=1-pchisq(QR,lag)
   ###
   cat("Rank-based Test: ","\n")
   cat("Test statistic: ",QR," p-value: ",pv2,"\n")
   
}


"MarchTest" <- function(zt,lag=10){
   #### Testing the presence of Multivariate ARCH effect
   if(!is.matrix(zt))zt=as.matrix(zt)
   nT=dim(zt)[1]; k=dim(zt)[2]
   C0=cov(zt)
   zt1=scale(zt,center=TRUE,scale=FALSE)  # subtract sample means
   C0iv=solve(C0)
   wk=zt1%*%C0iv
   wk=wk*zt1
   rt2=apply(wk,1,sum)-k
   m1=acf(rt2,lag.max=lag,plot=F)
   acf=m1$acf[2:(lag+1)]
   c1=c(1:lag)
   deno=rep(nT,lag)-c1
   Q=sum(acf^2/deno)*nT*(nT+2)
   pv1=1-pchisq(Q,lag)
   ###
   cat("Q(m) of squared series(LM test): ","\n")
   cat("Test statistic: ",Q," p-value: ",pv1,"\n")
   ##
   rk=rank(rt2)
   m2=acf(rk,lag.max=lag,plot=F)
   acf=m2$acf[2:(lag+1)]
   mu=-(rep(nT,lag)-c(1:lag))/(nT*(nT-1))
   v1=rep(5*nT^4,lag)-(5*c(1:lag)+9)*nT^3+9*(c(1:lag)-2)*nT^2+2*c(1:lag)*(5*c(1:lag)+8)*nT+16*c(1:lag)^2
   v1=v1/(5*(nT-1)^2*nT^2*(nT+1))
   QR=sum((acf-mu)^2/v1)
   pv2=1-pchisq(QR,lag)
   ###
   cat("Rank-based Test: ","\n")
   cat("Test statistic: ",QR," p-value: ",pv2,"\n")
   
   ### mq statistics
   cat("Q_k(m) of squared series: ","\n")
   x=zt^2
   g0=var(x)
   ginv=solve(g0)
   qm=0.0
   df = 0
   for (i in 1:lag){
      x1=x[(i+1):nT,]
      x2=x[1:(nT-i),]
      g = cov(x1,x2)
      g = g*(nT-i-1)/(nT-1)
      h=t(g)%*%ginv%*%g%*%ginv
      qm=qm+nT*nT*sum(diag(h))/(nT-i)
      df=df+k*k
   }
   pv3=1-pchisq(qm,df)
   cat("Test statistic: ",qm," p-value: ",pv3,"\n")
   #### Robust approach via trimming
   cut1=quantile(rt2,0.95)
   idx=c(1:nT)[rt2 <= cut1]
   x=zt[idx,]^2
   eT=length(idx)
   g0=var(x)
   ginv=solve(g0)
   qm=0.0
   df = 0
   for (i in 1:lag){
      x1=x[(i+1):eT,]
      x2=x[1:(eT-i),]
      g = cov(x1,x2)
      g = g*(eT-i-1)/(eT-1)
      h=t(g)%*%ginv%*%g%*%ginv
      qm=qm+eT*eT*sum(diag(h))/(eT-i)
      df=df+k*k
   }
   pv4=1-pchisq(qm,df)
   cat("Robust Test(5%) : ",qm," p-value: ",pv4,"\n")
   
}

############################
"BEKK11" <- function(rt,include.mean=T,cond.dist="normal"){
   # Estimation of a bivariate or 3-dimensional BEKK(1,1) model
   # rt: return series
   if(!is.matrix(rt))rt=as.matrix(rt)
   nT=dim(rt)[1]
   k=dim(rt)[2]
   if(k > 3){
      cat("Program Note: Dimension is limited to 3","\n")
      k=3
   }
   RTN <- rt[,1:k]
   #
   # obtain some initial estimates of the parameters
   mu=apply(RTN,2,mean)
   Cov1 <- cov(RTN); S=0.000001; S1=-0.5
   if(k==2){
      A1=matrix(c(0.1,0.02,0.02,0.1),k,k)
      B1=matrix(c(0.8,0.1,0.1,0.8),k,k)
      m1=chol(Cov1)
      #
      if(include.mean){
         par=c(mu1=mu[1],mu2=mu[2],A011=m1[1,1],A021=m1[1,2],A022=m1[2,2],A11=A1[1,1],A21=A1[2,1],A12=A1[1,2],A22=A1[2,2],
         B11=B1[1,1],B21=B1[2,1],B12=B1[1,2],B22=B1[2,2])
         c1=c(mu1=-10*abs(mu[1]),mu2=-10*abs(mu[2]),A011=m1[1,1]*0.2,A021=m1[1,2]*0.2,A022=m1[2,2]*0.2,A11=S,A21=S1,A12=S1,
         A22=S,B11=S,B21=S1,B12=S1,B22=S)
         c2=c(mu1=10*abs(mu[1]),mu2=10*abs(mu[2]),A011=m1[1,1]*1.1,A021=m1[1,2]*1.1,A022=m1[2,2]*1.1,A11=1-S,A21=-S1,A12=-S1,
         A22=1-S,B11=1-S,B21=-S1,B12=-S1,B22=1-S)
      }
      else{
         par=c(A011=m1[1,1],A021=m1[1,2],A022=m1[2,2],A11=A1[1,1],A21=A1[2,1],A12=A1[1,2],A22=A1[2,2],
         B11=B1[1,1],B21=B1[2,1],B12=B1[1,2],B22=B1[2,2])
         c1=c(A011=m1[1,1]*0.2,A021=m1[1,2]*0.2,A022=m1[2,2]*0.2,A11=S,A21=S1,A12=S1,
         A22=S,B11=S,B21=S1,B12=S1,B22=S)
         c2=c(A011=m1[1,1]*1.1,A021=m1[1,2]*1.1,A022=m1[2,2]*1.1,A11=1-S,A21=-S1,A12=-S1,
         A22=1-S,B11=1-S,B21=-S1,B12=-S1,B22=1-S)
      }
      ### end of if(k==2)
   }
   ##
   if(k==3){
      A1=matrix(c(0.1,0.02,0.02,0.02,0.1,0.02,0.02,0.02,0.1),k,k)
      B1=matrix(c(0.8,0.02,0.02,0.02,.8,0.02,0.02,0.02,0.8),k,k)
      m1=chol(Cov1)
      #
      if(include.mean){
         par=c(mu1=mu[1],mu2=mu[2],mu3=mu[3],A011=m1[1,1],A021=m1[1,2],A031=m1[1,3],A022=m1[2,2],A032=m1[2,3],A033=m1[3,3],A11=A1[1,1],A21=A1[2,1],A31=A1[3,1],A12=A1[1,2],A22=A1[2,2],A32=A1[3,2],A13=A1[1,3],A23=A1[2,3],A33=A1[3,3],
         B11=B1[1,1],B21=B1[2,1],B31=B1[3,1],B12=B1[1,2],B22=B1[2,2],B32=B1[3,2],B13=B1[1,3],B23=B1[2,3],B33=B1[3,3])
         c1=c(mu1=-10*abs(mu[1]),mu2=-10*abs(mu[2]),mu3=-10*abs(mu[3]),A011=m1[1,1]*0.2,A021=m1[1,2]*0.2,A031=m1[1,3]*0.2,A022=m1[2,2]*0.2,A032=m1[2,3]*0.2,A033=m1[3,3]*0.2,A11=S,A21=S1,A31=S1,A12=S1,
         A22=S,A32=S1,A13=S1,A23=S1,A33=S,B11=S,B21=S1,B31=S1,B12=S1,B22=S,B32=S1,B13=S1,B23=S1,B33=S)
         c2=c(mu1=10*abs(mu[1]),mu2=10*abs(mu[2]),mu3=10*abs(mu[3]),A011=m1[1,1]*1.1,A021=m1[1,2]*1.1,A031=m1[1,3]*1.1,A022=m1[2,2]*1.1,A032=m1[2,3]*1.1,A033=m1[3,3]*1.1,A11=1-S,A21=-S1,A31=-S1,A12=-S1,
         A22=1-S,A32=-S1,A13=-S1,A23=-S1,A33=1-S,B11=1-S,B21=-S1,B31=-S1,B12=-S1,B22=1-S,B32=-S1,B13=-S1,B23=-S1,B33=1-S)
      }
      else{
         par=c(A011=m1[1,1],A021=m1[1,2],A031=m1[1,3],A022=m1[2,2],A032=m1[2,3],A033=m1[3,3],A11=A1[1,1],A21=A1[2,1],A31=A1[3,1],A12=A1[1,2],A22=A1[2,2],A32=A1[3,2],A13=A1[1,3],A23=A1[2,3],A33=A1[3,3],
         B11=B1[1,1],B21=B1[2,1],B31=B1[3,1],B12=B1[1,2],B22=B1[2,2],B32=B1[3,2],B13=B1[1,3],B23=B1[2,3],B33=B1[3,3])
         c1=c(A011=m1[1,1]*0.2,A021=m1[1,2]*0.2,A031=m1[1,3]*0.2,A022=m1[2,2]*0.2,A032=m1[2,3]*0.2,A033=m1[3,3]*0.2,A11=S,A21=S1,A31=S1,A12=S1,
         A22=S,A32=S1,A13=S1,A23=S1,A33=S,B11=S,B21=S1,B31=S1,B12=S1,B22=S,B32=S1,B13=S1,B23=S1,B33=S)
         c2=c(A011=m1[1,1]*1.1,A021=m1[1,2]*1.1,A031=m1[1,3]*1.1,A022=m1[2,2]*1.1,A032=m1[2,3]*1.1,A033=m1[3,3]*1.1,A11=1-S,A21=-S1,A31=-S1,A12=-S1,
         A22=1-S,A32=-S1,A13=-S1,A23=-S1,A33=1-S,B11=1-S,B21=-S1,B31=-S1,B12=-S1,B22=1-S,B32=-S1,B13=-S1,B23=-S1,B33=1-S)
      }
      ### end of if(k==3)
   }
   cat("Initial estimates: ",par,"\n")
   cat("Lower limits: ",c1,"\n")
   cat("Upper limits: ",c2,"\n")
   #
   A1at <- function(x){
      # for multivariate GARCH models, this program compute the terms A1*at*t(at)*t(A1).
      x=as.matrix(x)
      Prod1=NULL
      for (i in 1:nrow(x)){
         Prod1=rbind(Prod1,kronecker(x[i,],x[i,]))
      }
      Prod1
   }
   #
   mlikeG <- function(par,RTN=RTN,include.mean=include.mean){
      # evaluate the log-likelihood function of a bivariate BEKK(1,1) model 
      # par: parameter vectors
      nT=dim(RTN)[1]; k=dim(RTN)[2]; Cov1=cov(RTN)
      #
      if(k==2){
         if(include.mean){
            mu1=par[1];mu2=par[2];A011=par[3];A021=par[4];A022=par[5];A11=par[6];A21=par[7];A12=par[8]
            A22=par[9];B11=par[10];B21=par[11];B12=par[12];B22=par[13]
         }
         else{
            mu1=0; mu2=0
            A011=par[1];A021=par[2];A022=par[3];A11=par[4];A21=par[5];A12=par[6]
            A22=par[7];B11=par[8];B21=par[9];B12=par[10];B22=par[11]
         }
         A0=matrix(c(A011,A021,0,A022),2,2)
         A0A0t=A0%*%t(A0)
         A1=matrix(c(A11,A21,A12,A22),2,2)
         B1=matrix(c(B11,B21,B12,B22),2,2)
         resi=cbind(RTN[,1]-mu1,RTN[,2]-mu2)
         Sig=matrix(c(Cov1),1,4)
         res=resi%*%t(A1)
         ArchP=A1at(res)
      }
      #
      if(k==3){
         if(include.mean){
            mu1=par[1];mu2=par[2];mu3=par[3];A011=par[4];A021=par[5];A031=par[6];A022=par[7];A032=par[8];A033=par[9];A11=par[10];A21=par[11];A31=par[12];A12=par[13];A22=par[14];A32=par[15];A13=par[16];A23=par[17];A33=par[18];B11=par[19];B21=par[20];B31=par[21];
            B12=par[22];B22=par[23];B32=par[24];B13=par[25];B23=par[26];B33=par[27]
         }
         else{
            mu1=0; mu2=0; mu3=0;
            A011=par[1];A021=par[2];A031=par[3];A022=par[4];A032=par[5];A033=par[6];A11=par[7];A21=par[8];A31=par[9];A12=par[10];A22=par[11];A32=par[12];A13=par[13];A23=par[14];A33=par[15];B11=par[16];B21=par[17];B31=par[18];
            B12=par[19];B22=par[20];B32=par[21];B13=par[22];B23=par[23];B33=par[24]
         }
         A0=matrix(c(A011,A021,A031,0,A022,A032,0,0,A033),k,k)
         A0A0t=A0%*%t(A0)
         A1=matrix(c(A11,A21,A31,A12,A22,A32,A13,A23,A33),k,k)
         B1=matrix(c(B11,B21,B31,B12,B22,B32,B13,B23,B33),k,k)
         resi=cbind(RTN[,1]-mu1,RTN[,2]-mu2,RTN[,3]-mu3)
         Sig=matrix(c(Cov1),1,k^2)
         res=resi%*%t(A1)
         ArchP=A1at(res)
      }
      #
      llike=0
      for (t in 2:nT){
         Sigt=A0A0t+matrix(ArchP[t-1,],k,k)+B1%*%matrix(Sig[t-1,],k,k)%*%t(B1)
         Sigt=(Sigt+t(Sigt))/2
         Sig=rbind(Sig,c(Sigt))
         d1=dmvnorm(resi[t,],mean=rep(0,k),Sigt,log=TRUE)
         llike=llike-d1
      }
      llike
   }
   # Estimate Parameters and Compute Numerically Hessian:
   fit = nlminb(start = par, objective = mlikeG, RTN=RTN, include.mean=include.mean, 
   lower = c1, upper = c2) ##, control = list(trace=3))
   epsilon = 0.0003 * fit$par
   if(k==3)epsilon=0.0005*fit$par
   npar=length(par)
   Hessian = matrix(0, ncol = npar, nrow = npar)
   for (i in 1:npar) {
      for (j in 1:npar) {
         x1 = x2 = x3 = x4  = fit$par
         x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
         x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
         x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
         x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
         Hessian[i, j] = (mlikeG(x1,RTN=RTN,include.mean=include.mean)-mlikeG(x2,RTN=RTN,include.mean=include.mean)
         -mlikeG(x3,RTN=RTN,include.mean=include.mean)+mlikeG(x4,RTN=RTN,include.mean=include.mean))/
         (4*epsilon[i]*epsilon[j])
      }
   }
   # Step 6: Create and Print Summary Report:
   est=fit$par
   se.coef = sqrt(diag(solve(Hessian)))
   tval = fit$par/se.coef
   matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
   dimnames(matcoef) = list(names(tval), c(" Estimate",
   " Std. Error", " t value", "Pr(>|t|)"))
   cat("\nCoefficient(s):\n")
   printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
   #
   BEKK11vol <- function(rtn,est,include.mean=T){
      ## Compute the volatility matrices of a fitted BEKK11 model.
      if(!is.matrix(rtn))rtn=as.matrix(rtn)
      nT=dim(rtn)[1]; k=dim(rtn)[2]
      if(k > 3){k=3; rtn=rtn[,1:k]}
      mu1=0; mu2=0; mu3=0; icnt=0
      if(k==2){
         if(include.mean){
            mu1=est[1]; mu2=est[2]; icnt=k
         }
         A0=matrix(c(est[icnt+1],est[icnt+2],0,est[icnt+3]),2,2)
         A0A0t=A0%*%A0
         A1=matrix(c(est[(icnt+4):(icnt+7)]),2,2)
         B1=matrix(c(est[(icnt+8):(icnt+11)]),2,2)
         resi=cbind(rtn[,1]-mu1,rtn[,2]-mu2)
         Cov1=cov(rtn)
         Sig=matrix(c(Cov1),1,4)
         res=resi%*%t(A1)
         ArchP=A1at(res)
      }
      if(k==3){
         if(include.mean){
            mu1=est[1]; mu2=est[2]; mu3=est[3]; icnt=k
         }
         A0=matrix(c(est[icnt+1],est[icnt+2],est[icnt+3],0,est[icnt+4],est[icnt+5],0,0,est[icnt+6]),k,k)
         A0A0t=A0%*%t(A0)
         icnt=icnt+6
         A1=matrix(c(est[(icnt+1):(icnt+9)]),k,k)
         icnt=icnt+9
         B1=matrix(c(est[(icnt+1):(icnt+9)]),k,k)
         resi=cbind(rtn[,1]-mu1,rtn[,2]-mu2,rtn[,3]-mu3)
         Cov1=cov(rtn)
         Sig=matrix(c(Cov1),1,k^2)
         res=resi%*%t(A1)
         ArchP=A1at(res)
      }
      #
      for (t in 2:nT){
         Sigt=A0A0t+matrix(ArchP[t-1,],k,k)+B1%*%matrix(Sig[t-1,],k,k)%*%t(B1)
         Sigt=(Sigt+t(Sigt))/2
         Sig=rbind(Sig,c(Sigt))
      }
      BEKK11vol <- list(volmtx=Sig)
   }
   
   m2=BEKK11vol(RTN,est,include.mean=include.mean)
   Sigma.t=m2$volmtx
   
   BEKK11 <- list(estimates=est,HessianMtx=Hessian,Sigma.t=Sigma.t)
}



