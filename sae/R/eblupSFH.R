eblupSFH <-
function(formula,vardir,proxmat,method="REML",MAXITER=100,PRECISION=0.0001,data)
{
   result <- list(eblup=NA, 
                  fit=list(method=method, convergence=TRUE, iterations=0, estcoef=NA, 
                  refvar=NA, spatialcorr=NA, goodness=NA)
                 ) 

   if (method!="REML" & method!="ML")
       stop(" method=\"",method, "\" must be \"REML\" or \"ML\".")

   namevar     <- deparse(substitute(vardir))
   if (!missing(data))
   {
      formuladata <- model.frame(formula,na.action = na.omit,data)
      X           <- model.matrix(formula,data)        
      vardir      <- data[,namevar]
   } else
   {
      formuladata <- model.frame(formula,na.action = na.omit)
      X           <- model.matrix(formula)        
   }
   y <- formuladata[,1]            

   if (attr(attributes(formuladata)$terms,"response")==1)
      textformula <- paste(formula[2],formula[1],formula[3])
   else
      textformula <- paste(formula[1],formula[2])

   if (length(na.action(formuladata))>0)
      stop("Argument formula=",textformula," contains NA values.")
   if (any(is.na(vardir)))
      stop("Argument vardir=",namevar," contains NA values.")

   proxmatname <- deparse(substitute(proxmat))
   if (any(is.na(proxmat)))
      stop("Argument proxmat=",proxmatname," contains NA values.")

   if (!is.matrix(proxmat))
      proxmat <- as.matrix(proxmat)

   nformula  <- nrow(X)   
   nvardir   <- length(vardir) 
   nproxmat  <- nrow(proxmat)
   if (nformula!=nvardir | nformula!=nproxmat)
      stop("   formula=",textformula," [rows=",nformula,"],\n", 
           "     vardir=",namevar," [rows=",nvardir,"] and \n",
           "     proxmat=",proxmatname," [rows=",nproxmat,"]\n",
           "  must be the same length.")
   if (nproxmat!=ncol(proxmat))
      stop(" Argument proxmat=",proxmatname," is not a square matrix [rows=",nproxmat,",columns=",ncol(proxmat),"].")

   m <-length(y)  # Sample size or number of areas
   p <-dim(X)[2]  # Num. of X columns of num. of auxiliary variables (including intercept)
   Xt <-t(X)
   yt <-t(y)
   proxmatt<-t(proxmat)
   I <-diag(1,m)
  
   # Initialize vectors containing estimators of variance and spatial correlation
   par.stim <-matrix(0,2,1)
   stime.fin<-matrix(0,2,1)
  
   # Initialize scores vector and Fisher information matrix
   s<-matrix(0,2,1)
   Idev<-matrix(0,2,2)

   # Initial value of variance set to the mean of sampling variances vardir
   # Initial value of spatial correlation set to 0.5
   sigma2.u.stim.S<-0
   rho.stim.S<-0
  
   sigma2.u.stim.S[1]<-median(vardir)
   rho.stim.S[1]<-0.5
  
   if (method=="REML")   # Fisher-scoring algorithm for REML estimators start
   { 
      k<-0
      diff.S<-PRECISION+1
      while ((diff.S>PRECISION)&(k<MAXITER))
      {
         k<-k+1

         # Derivative of covariance matrix V with respect to variance
         derSigma<-solve((I-rho.stim.S[k]*proxmatt)%*%(I-rho.stim.S[k]*proxmat))

         # Derivative of covariance matrix V with respect to spatial correlation parameter
         derRho<-2*rho.stim.S[k]*proxmatt%*%proxmat-proxmat-proxmatt
         derVRho<-(-1)*sigma2.u.stim.S[k]*(derSigma%*%derRho%*%derSigma)

         # Covariance matrix and inverse covariance matrix
         V<-sigma2.u.stim.S[k]*derSigma+I*vardir
         Vi<-solve(V)

         # Matrix P and coefficients'estimator beta
         XtVi<-Xt%*%Vi
         Q<-solve(XtVi%*%X)
         P<-Vi-t(XtVi)%*%Q%*%XtVi
         b.s<-Q%*%XtVi%*%y

         # Terms involved in scores vector and Fisher information matrix
         PD<-P%*%derSigma
         PR<-P%*%derVRho
         Pdir<-P%*%y

         # Scores vector
         s[1,1]<-(-0.5)*sum(diag(PD))+(0.5)*(yt%*%PD%*%Pdir)
         s[2,1]<-(-0.5)*sum(diag(PR))+(0.5)*(yt%*%PR%*%Pdir)

         # Fisher information matrix
         Idev[1,1]<-(0.5)*sum(diag(PD%*%PD))
         Idev[1,2]<-(0.5)*sum(diag(PD%*%PR))
         Idev[2,1]<-Idev[1,2]
         Idev[2,2]<-(0.5)*sum(diag(PR%*%PR))
      
         # Updating equation
         par.stim[1,1]<-sigma2.u.stim.S[k]
         par.stim[2,1]<-rho.stim.S[k]
      
         stime.fin<-par.stim+solve(Idev)%*%s

         # Restricting the spatial correlation to (-0.999,0.999)
         if (stime.fin[2,1]<=-1)
            stime.fin[2,1] <- -0.999
         if (stime.fin[2,1]>=1)
            stime.fin[2,1] <- 0.999
      
         # Restricting the spatial correlation to (-0.999,0.999) and the variance to (0.0001,infty)
         #if ((stime.fin[2,1]<=0.999)&(stime.fin[2,1]>=-0.999)&(stime.fin[1,1]>0.0001)){
            sigma2.u.stim.S[k+1]<-stime.fin[1,1]
            rho.stim.S[k+1]<-stime.fin[2,1]
            diff.S<-max(abs(stime.fin-par.stim)/par.stim)
         #}else
         #{
         #   sigma2.u.stim.S[k+1]<-stime.fin[1,1]
         #    rho.stim.S[k+1]<-stime.fin[2,1]
         #    diff.S<-PRECISION/10
         #}
      } # End of while
      
   } else       # Fisher-scoring algorithm for ML estimators start
   {
       k<-0
       diff.S<-PRECISION+1
       while ((diff.S>PRECISION)&(k<MAXITER))
       {
          k<-k+1

          # Derivative of covariance matrix V with respect to variance
          derSigma<-solve((I-rho.stim.S[k]*proxmatt)%*%(I-rho.stim.S[k]*proxmat))

          # Derivative of covariance matrix V with respect to spatial correlation
          derRho<-2*rho.stim.S[k]*proxmatt%*%proxmat-proxmat-proxmatt
          derVRho<-(-1)*sigma2.u.stim.S[k]*(derSigma%*%derRho%*%derSigma)

          # Covariance matrix and inverse covariance matrix
          V<-sigma2.u.stim.S[k]*derSigma+I*vardir
          Vi<-solve(V)

          # Coefficients'estimator beta and matrix P
          XtVi<-Xt%*%Vi
          Q<-solve(XtVi%*%X)
          P<-Vi-t(XtVi)%*%Q%*%XtVi
          b.s<-Q%*%XtVi%*%y

          # Terms involved in scores vector and Fisher information matrix
          PD<-P%*%derSigma
          PR<-P%*%derVRho
          Pdir<-P%*%y
          ViD<-Vi%*%derSigma
          ViR<-Vi%*%derVRho

          # Scores vector
          s[1,1]<-(-0.5)*sum(diag(ViD))+(0.5)*(yt%*%PD%*%Pdir)
          s[2,1]<-(-0.5)*sum(diag(ViR))+(0.5)*(yt%*%PR%*%Pdir)

          # Fisher information matrix
          Idev[1,1]<-(0.5)*sum(diag(ViD%*%ViD))
          Idev[1,2]<-(0.5)*sum(diag(ViD%*%ViR))
          Idev[2,1]<-Idev[1,2]
          Idev[2,2]<-(0.5)*sum(diag(ViR%*%ViR))
      
          # Updating equation      
          par.stim[1,1]<-sigma2.u.stim.S[k]
          par.stim[2,1]<-rho.stim.S[k]
      
          stime.fin<-par.stim+solve(Idev)%*%s

          # Restricting the spatial correlation to (-0.999,0.999)
          if (stime.fin[2,1]<=-1)
             stime.fin[2,1] <- -0.999
          if (stime.fin[2,1]>=1)
             stime.fin[2,1] <- 0.999
      
          # Restricting the spatial correlation to (-0.999,0.999) and the variance to (0.0001,infty)
          #if ((stime.fin[2,1]<=0.999)&(stime.fin[2,1]>=-0.999)&(stime.fin[1,1]>0.0001)){
             sigma2.u.stim.S[k+1]<-stime.fin[1,1]
             rho.stim.S[k+1]<-stime.fin[2,1]
             diff.S<-max(abs(stime.fin-par.stim)/par.stim)
          #}else
          #{
          #   sigma2.u.stim.S[k+1]<-stime.fin[1,1]
          #   rho.stim.S[k+1]<-stime.fin[2,1]
          #   diff.S<-PRECISION/10
          #}
    
      } # End of while

   }   
   
   # Final values of estimators
   if (rho.stim.S[k+1]==-0.999)
      rho.stim.S[k+1] <- -1
   else if (rho.stim.S[k+1]==0.999)
      rho.stim.S[k+1] <- 1
   rho <-rho.stim.S[k+1]

   sigma2.u.stim.S[k+1]<-max(sigma2.u.stim.S[k+1],0)
   sigma2u <- sigma2.u.stim.S[k+1]

   #print(rho.stim.S)
   #print(sigma2.u.stim.S)
           
   # Indicator of convergence
   result$fit$iterations  <- k  
   if(k>=MAXITER && diff>=PRECISION) 
   {
      result$fit$convergence <- FALSE
      return(result)
   }

   result$fit$refvar       <- sigma2u
   result$fit$spatialcorr <- rho
   if (sigma2u<0 || rho<(-1) || rho>1 )  # COMPROBAR
   {
      print("eblupSFH: este mensaje no debe salir")
      return(result)
   }
   # Computation of the coefficients'estimator (Bstim)
   A   <-solve((I-rho*proxmatt)%*%(I-rho*proxmat))    
   G   <-sigma2u*A
   V   <-G+I*vardir
   Vi  <-solve(V)
   XtVi<-Xt%*%Vi
   Q   <-solve(XtVi%*%X)
   Bstim<-Q%*%XtVi%*%y
 
   # Significance of the regression coefficients
   std.errorbeta<-sqrt(diag(Q))
   tvalue       <-Bstim/std.errorbeta
   pvalue       <-2*pnorm(abs(tvalue),lower.tail=FALSE)
   coef         <-data.frame(beta=Bstim,std.error=std.errorbeta,tvalue,pvalue)
  
   # Goodness of fit measures: loglikelihood, AIC, BIC
   Xbeta <-X%*%Bstim
   resid <-y-Xbeta
   loglike<-(-0.5)*(m*log(2*pi)+determinant(V,logarithm=TRUE)$modulus+t(resid)%*%Vi%*%resid)
   AIC    <-(-2)*loglike+2*(p+2)
   BIC    <-(-2)*loglike+(p+2)*log(m)
   goodness<-c(loglike=loglike,AIC=AIC,BIC=BIC)

   # Computation of the Spatial EBLUP
   res<-y-X%*%Bstim
   thetaSpat<-X%*%Bstim+G%*%Vi%*%res

   result$fit$estcoef  <- coef
   result$fit$goodness <- goodness
   result$eblup       <- thetaSpat

   return(result)
}