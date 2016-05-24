npbmseSFH <-
function(formula,vardir,proxmat,B=100,method="REML",MAXITER=100,PRECISION=0.0001,data)
{   
   result <- list(est=NA, mse=NA)

   if (method!="REML")
      stop(" method=\"",method, "\" must be \"REML\".")

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
      stop("Argument proxmat=",proxmatname," is not a square matrix [rows=",nproxmat,",columns=",ncol(proxmat),"].")



   m<-dim(X)[1]  # Sample size or number of areas
   p<-dim(X)[2]  # Num. of auxiliary variables (including intercept)

   # Fit the model to initial sample data using the given method
   result$est<-try(eblupSFH(y~X-1,vardir,proxmat,method,MAXITER,PRECISION)) 
   if (result$est$fit$convergence==FALSE) 
   {
      warning("The fitting method does not converge.\n") 
      return (result);
   }

   # Initial estimators of model coefficients, variance and spatial 
   # correlation which will act as true values in the bootstrap 
   # procedure.

   Bstim.boot <-result$est$fit$estcoef$beta          
   rho.boot   <-result$est$fit$spatialcorr
   sigma2.boot<-result$est$fit$refvar           

   # Auxiliary calculations

   I<-diag(1,m)
   proxmatt<-t(proxmat)
   Xt<-t(X)

   Irhoproxmat<-I-rho.boot*proxmat
   Irhoproxmatt<-t(Irhoproxmat)
   Ar<-solve(Irhoproxmatt%*%Irhoproxmat)
   Gr<-sigma2.boot*Ar
   Vr<-Gr+I*vardir
   Vri<-solve(Vr)
   Qr<-solve(Xt%*%Vri%*%X)

   # Analytical estimators of g1 and g2, used for the bias-corrected 
   # PB MSE estimator.

   g1sp<-rep(0,m)
   g2sp<-rep(0,m)

   XtVri<-Xt%*%Vri
   Qr<-solve(XtVri%*%X)

   # Calculate g1 and g2

   Ga<-Gr-Gr%*%Vri%*%Gr
   Gb<-Gr%*%t(XtVri)
   Xa<-matrix(0,1,p)

   for (i in 1:m) {
     g1sp[i]<-Ga[i,i]
     Xa[1,]<-X[i,]-Gb[i,]
     g2sp[i]<-Xa%*%Qr%*%t(Xa)
   }

   # Residual vectors
   res<-y-X%*%Bstim.boot
   vstim<-Gr%*%Vri%*%res

   # Calculate covariance matrices of residual vectors
   VG<-Vr-Gr
   P<-Vri-Vri%*%X%*%Qr%*%Xt%*%Vri

   Ve<-VG%*%P%*%VG
   Vu<-Irhoproxmat%*%Gr%*%P%*%Gr%*%Irhoproxmatt

   # Square roots of covariance matrices

   VecVe0<-eigen(Ve)$vectors
   VecVe<-VecVe0[,1:(m-p)]
   ValVe0<-eigen(Ve)$values
   Valve<-diag(sqrt(1/ValVe0[1:(m-p)]))
   Vei05<-VecVe%*%Valve%*%t(VecVe)

   VecVu0<-eigen(Vu)$vectors
   VecVu<-VecVu0[,1:(m-p)]
   ValVu0<-1/(eigen(Vu)$values)
   ValVu<-diag(sqrt(ValVu0[1:(m-p)]))
   Vui05<-VecVu%*%ValVu%*%t(VecVu)

   # Standardize residual vectors
   ustim<-as.vector(Vui05%*%((Irhoproxmat)%*%vstim))
   estim<-as.vector(Vei05%*%(res-vstim))
   sdu<-sqrt(sigma2.boot)

   u.std<-rep(0,m)
   e.std<-rep(0,m)

   for (i in 1:m){
     u.std[i]<-(sdu*(ustim[i]-mean(ustim)))/
     sqrt(mean((ustim-mean(ustim))^2))
     e.std[i]<-(estim[i]-mean(estim))/
     sqrt(mean((estim-mean(estim))^2))
   }

   # Bootstrap algorithm starts

   difmse.npb<-matrix(0,m,1)
   difg3Spat.npb<-matrix(0,m,1)
   g1sp.aux<-matrix(0,m,1)
   g2sp.aux<-matrix(0,m,1)
   difg1sp.npb<-matrix(0,m,1)
   difg2sp.npb<-matrix(0,m,1)


   cat("\nBootstrap procedure with B =",B,"iterations starts.\n")
   boot <- 1
   while (boot<=B) 
   {
      # Generate boostrap data
      u.boot     <-sample(u.std,m,replace=TRUE)
      e.samp     <-sample(e.std,m,replace=TRUE)
      e.boot     <-sqrt(vardir)*e.samp
      v.boot     <-solve(Irhoproxmat)%*%u.boot
      theta.boot <-X%*%Bstim.boot+v.boot
      direct.boot<-theta.boot+e.boot

      # Fit the model to bootstrap data
      results.SpFH.boot<-eblupSFH(direct.boot[,1]~X-1,vardir,proxmat,method,MAXITER,PRECISION)

      # Generate a new sample if estimators are not satisfactory
      if (results.SpFH.boot$fit$convergence==FALSE | results.SpFH.boot$fit$refvar<0 | 
          results.SpFH.boot$fit$spatialcorr<(-1) | results.SpFH.boot$fit$spatialcorr>1)
         next

      cat("b =",boot,"\n")

      Bstim.ML.boot<-results.SpFH.boot$fit$estcoef[,1]                                   
      rho.ML.boot<-results.SpFH.boot$fit$spatialcorr
      sigma2.ML.boot<-results.SpFH.boot$fit$refvar                                        
      thetaEB.SpFH.boot<-results.SpFH.boot$eblup                                   

      # Nonparametric bootstrap estimator of g3
      Bstim.sblup<-Qr%*%XtVri%*%direct.boot[,1]
      thetaEB.SpFH.sblup.boot<-X%*%Bstim.sblup+Gr%*%Vri%*%
      (direct.boot[,1]-X%*%Bstim.sblup)

      difg3Spat.npb[,1]<-difg3Spat.npb[,1]+
      (thetaEB.SpFH.boot-thetaEB.SpFH.sblup.boot)^2

      # Naive nonparametric bootstrap MSE
      difmse.npb[,1]<-difmse.npb[,1]+(thetaEB.SpFH.boot[,1]-theta.boot)^2

      # g1 and g2 for each bootstrap sample
      A<-solve((I-rho.ML.boot*proxmatt)%*%(I-rho.ML.boot*proxmat))
      G<-sigma2.ML.boot*A
      V<-G+I*vardir
      Vi<-solve(V)
      XtVi<-Xt%*%Vi
      Q<-solve(XtVi%*%X)

      Ga<-G-G%*%Vi%*%G
      Gb<-G%*%Vi%*%X
      Xa<-matrix(0,1,p)

      for (i in 1:m){
         g1sp.aux[i]<-Ga[i,i] 
         Xa[1,]<-X[i,]-Gb[i,]
         g2sp.aux[i]<-Xa%*%Q%*%t(Xa)
      }

      difg1sp.npb<-difg1sp.npb+g1sp.aux
      difg2sp.npb<-difg2sp.npb+g2sp.aux

      boot <- boot+1
   } # End of bootstrap cycle

   # Final naive nonparametric bootstrap MSE estimator
   mse.npb<-difmse.npb[,1]/B

   # Final bias-corrected nonparametric bootstrap MSE estimator
   g3Spat.npb<-difg3Spat.npb/B
   g1sp.npb<-difg1sp.npb/B
   g2sp.npb<-difg2sp.npb/B
   mse.npb2<-2*(g1sp+g2sp)-difg1sp.npb[,1]/B-difg2sp.npb[,1]/B+
   difg3Spat.npb[,1]/B

   result$mse   <- data.frame(mse=mse.npb, msebc=mse.npb2)
   return(result)
}

