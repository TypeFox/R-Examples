mseSFH <-
function(formula,vardir,proxmat,method="REML",MAXITER=100,PRECISION=0.0001,data)
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

   A    <- result$est$fit$refvar
   rho  <- result$est$fit$spatialcorr

   m<-dim(X)[1]  # Sample size or number of areas
   p<-dim(X)[2]  # Number of auxiliary variables

   # Initialize vectors containing the values of g1-g4 and mse for each area
   g1d<-rep(0,m)
   g2d<-rep(0,m)
   g3d<-rep(0,m)
   g4d<-rep(0,m)
   mse2d.aux<-rep(0,m)
   mse2d<-rep(0,m)
  
   # Auxiliary calculations
   I<-diag(1,m)
   Xt<-t(X)
   proxmatt<-t(proxmat)
   Ci<-solve((I-rho*proxmatt)%*%(I-rho*proxmat))
   G<-A*Ci
   V<-G+I*vardir
   Vi<-solve(V)
   XtVi<-Xt%*%Vi
   Q<-solve(XtVi%*%X)
  
   Ga<-G-G%*%Vi%*%G
  
   Gb<-G%*%Vi%*%X
   Xa<-matrix(0,1,p)
  
   for (i in 1:m) { 
      g1d[i]<-Ga[i,i]  # g1 contains the diagonal elements of Ga
      Xa[1,]<-X[i,]-Gb[i,]
      g2d[i]<-Xa%*%Q%*%t(Xa)
   }
  
   # Auxiliary calculations for g3
   derRho<-2*rho*proxmatt%*%proxmat-proxmat-proxmatt
   Amat<-(-1)*A*(Ci%*%derRho%*%Ci)
   P<-Vi-t(XtVi)%*%Q%*%XtVi
   PCi<-P%*%Ci
   PAmat<-P%*%Amat
  
   Idev<-matrix(0,2,2)
   Idev[1,1]<-(0.5)*sum(diag((PCi%*%PCi)))
   Idev[1,2]<-(0.5)*sum(diag((PCi%*%PAmat)))
   Idev[2,1]<-Idev[1,2]
   Idev[2,2]<-(0.5)*sum(diag((PAmat%*%PAmat)))
   Idevi<-solve(Idev)
  
   ViCi<-Vi%*%Ci
   ViAmat<-Vi%*%Amat
  
   # Calculation of g3
   l1<-ViCi-A*ViCi%*%ViCi
   l1t<-t(l1)
   l2<-ViAmat-A*ViAmat%*%ViCi
   l2t<-t(l2)
   L<-matrix(0,2,m)
   for (i in 1:m)
   {
      L[1,]<-l1t[i,]
      L[2,]<-l2t[i,]
      g3d[i]<-sum(diag(L%*%V%*%t(L)%*%Idevi))
   }
  
   mse2d.aux<-g1d+g2d+2*g3d
  
   # Bias correction of Singh et al
   psi<-diag(vardir,m)
   D12aux<-(-1)*(Ci%*%derRho%*%Ci)
   D22aux<-2*A*Ci%*%derRho%*%Ci%*%derRho%*%Ci-2*A*Ci%*%proxmatt%*%proxmat%*%Ci
   D<-(psi%*%Vi%*%D12aux%*%Vi%*%psi)*(Idevi[1,2]+Idevi[2,1])+psi%*%Vi%*%D22aux%*%Vi%*%psi*Idevi[2,2]
  
   for (i in 1:m) {g4d[i]<-(0.5)*D[i,i]}
  
   # Computation of estimated MSE of Singh et al
   mse2d<-mse2d.aux-g4d

   if (method=="ML")
   {
      # Calculate bML
      QXtVi <- Q%*%XtVi
      ViX   <- Vi%*%X
      h1    <- (-1)*sum(diag(QXtVi%*%Ci%*%ViX))
      h2    <- (-1)*sum(diag(QXtVi%*%Amat%*%ViX))
      h     <- matrix(c(h1,h2),nrow=2,ncol=1)
      bML   <- (Idevi%*%h)/2 
      tbML  <- t(bML)

      # Calculate gradient of g1d
      GVi     <- G%*%Vi   # G<-A*Ci
      GViCi   <- GVi%*%Ci
      GViAmat <- GVi%*%Amat
      ViCi    <- Vi%*%Ci
      dg1_dA  <- Ci   - 2*GViCi   + A*GViCi%*%ViCi
      dg1_dp  <- Amat - 2*GViAmat + A*GViAmat%*%ViCi
      gradg1d <- matrix(0,nrow=2,ncol=1) 

      bMLgradg1 <- rep(0,m)   
      for (i in 1:m)
      {
         gradg1d[1,1] <- dg1_dA[i,i]
         gradg1d[2,1] <- dg1_dp[i,i]
         bMLgradg1[i] <- tbML%*%gradg1d             
      }
      mse2d <- mse2d - bMLgradg1
   }

   result$mse <- mse2d   
   return (result)
}

