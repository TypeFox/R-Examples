pbmseSFH <-
function(formula,vardir,proxmat,B=100,method="REML",MAXITER=100,PRECISION=0.0001,data)
{ 
   result <- list(est=NA, mse=NA)

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


   result$est <- eblupSFH(y~X-1,vardir,proxmat,method,MAXITER,PRECISION)
   if (result$est$fit$convergence==FALSE) 
   {
      warning("The fitting method does not converge.\n") 
      return (result);
   }

   beta <- result$est$fit$estcoef$beta
   A    <- result$est$fit$refvar
   rho  <- result$est$fit$spatialcorr
  
   m<-dim(X)[1]  # Sample size or number of areas
   p<-dim(X)[2]  # Num. of auxiliary variables (including intercept)

   # Initial estimators of model coefficients, variance and spatial correlation
   # which will act as true values in the bootstrap procedure.
   Bstim.boot<-beta
   rho.boot  <-rho
   sigma2.boot<-A

   I  <-diag(1,m)
   Xt <-t(X)
   proxmatt <-t(proxmat)
  
   # Analytical estimators of g1 and g2, used for the bias-corrected PB MSE estimator
   g1sp<-rep(0,m)
   g2sp<-rep(0,m)
  
   Amat.sblup<-solve((I-rho.boot*proxmatt)%*%(I-rho.boot*proxmat))
   G.sblup<-sigma2.boot*Amat.sblup
   V.sblup<-G.sblup+I*vardir
   V.sblupi<-solve(V.sblup)
   XtV.sblupi<-Xt%*%V.sblupi
   Q.sblup<-solve(XtV.sblupi%*%X)
  
   # Calculate g1
   Ga.sblup<-G.sblup-G.sblup%*%V.sblupi%*%G.sblup
   for (i in 1:m) {g1sp[i]<-Ga.sblup[i,i]}

   # Calculate g2
   Gb.sblup<-G.sblup%*%t(XtV.sblupi)
   Xa.sblup<-matrix(0,1,p)
   for (i in 1:m){
      Xa.sblup[1,]<-X[i,]-Gb.sblup[i,]
      g2sp[i]<-Xa.sblup%*%Q.sblup%*%t(Xa.sblup)
   }

   # Initialize vectors adding g1, g2, g3 and naive PB MSE estimators for each bootstrap replication
   summse.pb<-rep(0,m)
   sumg1sp.pb<-rep(0,m)
   sumg2sp.pb<-rep(0,m)
   sumg3sp.pb<-rep(0,m)
   g1sp.aux<-rep(0,m)
   g2sp.aux<-rep(0,m)

   # Bootstrap cycle starts
   cat("\nBootstrap procedure with B =",B,"iterations starts.\n")
   boot <- 1
   while (boot<=B) 
   {
      # Generate a bootstrap sample
      u.boot     <-rnorm(m,0,sqrt(sigma2.boot))
      v.boot     <-solve(I-rho.boot*proxmat)%*%u.boot
      theta.boot <-X%*%Bstim.boot+v.boot
      e.boot     <-rnorm(m,0,sqrt(vardir))
      direct.boot<-theta.boot+e.boot

      # Fit of the model to bootstrap data
      resultsSp  <-eblupSFH(direct.boot~X-1,vardir,proxmat,method,MAXITER,PRECISION)

      # Generate a new sample if estimators are not satisfactory
      # As variance is >=0 and correlation is within the range, this message is not printed
      if (resultsSp$fit$convergence==FALSE | resultsSp$fit$refvar<0 | resultsSp$fit$spatialcorr<(-1) | resultsSp$fit$spatialcorr>1)
      {
         print("eblupSFH: If this message is printed, please contact with the authors of the package sae.")
         next
      }
      cat("b =",boot,"\n")

      # Fit of the model to bootstrap data
      sigma2.simula.ML<-resultsSp$fit$refvar                         
      rho.simula.ML   <-resultsSp$fit$spatialcorr
      beta.ML         <-resultsSp$fit$estcoef$beta                           

      # Calculation of the bootstrap Spatial EBLUP
      Amat<-solve((I-rho.simula.ML*proxmatt)%*%(I-rho.simula.ML*proxmat))
      G<-sigma2.simula.ML*Amat
      V<-G+I*vardir
      Vi<-solve(V)
      Xbeta<-X%*%beta.ML
      thetaEBLUPSpat.boot<-Xbeta+G%*%Vi%*%(direct.boot-Xbeta)

      # Naive parametric bootstrap MSE
      summse.pb<-summse.pb+(thetaEBLUPSpat.boot-theta.boot)^2

      # Bias-corrected parametric bootstrap:
      # For de bias of g1 and g2, calculate g1sp and g2sp for each bootstrap sample
      XtVi<-Xt%*%Vi
      Q<-solve(XtVi%*%X)

      Ga<-G-G%*%Vi%*%G
      for (i in 1:m) {g1sp.aux[i]<-Ga[i,i]} # g1sp contains the diagonal elements of Ga

      Gb<-G%*%Vi%*%X
      Xa<-matrix(0,1,p)
      for (i in 1:m){
         Xa[1,]<-X[i,]-Gb[i,]
         g2sp.aux[i]<-Xa%*%Q%*%t(Xa)
      }

      # Bootstrap spatial BLUP
      Bstim.sblup<-solve(XtV.sblupi%*%X)%*%XtV.sblupi%*%direct.boot
      Xbeta.sblup<-X%*%Bstim.sblup
      thetaEBLUPSpat.sblup.boot<-Xbeta.sblup+G.sblup%*%V.sblupi%*%(direct.boot-Xbeta.sblup)

      # Parametric bootstrap estimator of g3
      sumg3sp.pb<-sumg3sp.pb+(thetaEBLUPSpat.boot-thetaEBLUPSpat.sblup.boot)^2

      # Expectation of estimated g1 and g2
      sumg1sp.pb<-sumg1sp.pb+g1sp.aux
      sumg2sp.pb<-sumg2sp.pb+g2sp.aux

      boot <- boot+1
   } # End of bootstrap cycle

   # Final naive parametric bootstrap MSE estimator
   mse.pb<-summse.pb/B

   # Final bias-corrected bootstrap MSE estimator
   g1sp.pb<-sumg1sp.pb/B
   g2sp.pb<-sumg2sp.pb/B
   g3sp.pb<-sumg3sp.pb/B
   mse.pb2<-2*(g1sp+g2sp)-g1sp.pb-g2sp.pb+g3sp.pb

  # Return naive and bias-corrected parametric bootstrap
  result$mse   <- data.frame(mse=mse.pb, msebc=mse.pb2)

  return(result)
}

