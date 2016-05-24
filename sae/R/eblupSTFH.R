eblupSTFH <-
function(formula,D,T,vardir,proxmat,model="ST",MAXITER=100,PRECISION=0.0001,data)
{
   result <- list(eblup=NA, 
                  fit=list(model=model, convergence=TRUE, iterations=0, estcoef=NA, 
                  estvarcomp=NA, goodness=NA)
                 ) 

   if (model!="S" && model!="ST")
      stop("Argument model=\"",model,"\" must be \"S\" or \"ST\".")

   namevar <- deparse(substitute(vardir))
   if (!missing(data))
   {
      formuladata <- model.frame(formula,na.action = na.omit,data)
      X           <- model.matrix(formula,data)        
      vardir    <- data[,namevar]
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

   nformula  <- nrow(X)   
   nvardir   <- length(vardir) 
   M <- D*T
   if (nformula!=nvardir | nformula!=M)
      stop("   formula=",textformula," [rows=",nformula,"] and \n", 
           "     vardir=",namevar," [rows=",nvardir,"] \n",
           "  must have D*T = ",D,"*",T," = ",M," rows.")

   if (!is.matrix(proxmat))
      proxmat <- as.matrix(proxmat)
   nproxmat  <- nrow(proxmat)
   if (nproxmat!=D || ncol(proxmat)!=D)
      stop("Argument proxmat=",proxmatname," [rows=",nproxmat,",columns=",ncol(proxmat),
           "] must be a square matrix of size D=",D,".")

                                        # Initialize theta0 values
   sigma2_0 <- 0.5*median(vardir)
   if (model=="S")
      theta0 <- rbind(sigma21=sigma2_0, rho1=0.5, sigma22=sigma2_0)
   else
      theta0 <- rbind(sigma21=sigma2_0, rho1=0.5, sigma22=sigma2_0, rho2=0.5)
   nparam <- nrow(theta0)

                                      # Initialization
   invA   <- matrix(0,nrow=M,ncol=M)
   invAZ1 <- matrix(0,nrow=M,ncol=D)
   tZ1PZ1 <- matrix(0,nrow=D, ncol=D) 
   tZ1P   <- matrix(0,nrow=D,ncol=M)  
   S <- trPV   <- matrix(0,nrow=nparam,ncol=1)
   F <- trPVPV <- matrix(0,nrow=nparam,ncol=nparam)

   Va     <- list()
                                       ### Calculate Z1
   vector1T <- matrix(1,nrow=T, ncol=1)
   Z1       <- matrix(0,nrow=D*T, ncol=D)
   first <- 1
   for (d in 1:D)
   {
      last <- first+T-1
      Z1[first:last,d]<-vector1T
      first <- last+1
   } 

   tZ1 <- t(Z1)                   
   ty  <- t(y)
   tX  <- t(X)
   tproxmatproxmat <- crossprod(proxmat)
   p1derivrho1 <- - proxmat - t(proxmat) 
   Id      <- diag(1,nrow=D,ncol=D)
   Tmen1   <- T-1

   if (model=="S")
      Va[[3]] <- diag(1,nrow=M,ncol=M)                          
   else
   {
      PV <- list()                                                   
      Omega2drho2 <- derivOmega2drho2 <- matrix(0,nrow=T,ncol=T)     
      seqTmen1_1  <- sequence(Tmen1:1)                                
      Ve <- diag(vardir)                                           
   }

   thetakmas1 <- thetak <- theta0
   k <- 0
   diff <- PRECISION+1
   while (diff>PRECISION & k<MAXITER)
   {
      k <- k+1
      thetak    <- thetakmas1
      sigma21_k <- thetak["sigma21",1]
      rho1_k    <- thetak["rho1",1]
      sigma22_k <- thetak["sigma22",1]

      Omega1rho1_k <- try(solve(crossprod(Id-rho1_k*proxmat)))  
      if(class(Omega1rho1_k)=="try-error")
      {
         result$fit$convergence <- FALSE
         return (result)
      }
        
      Vu1 <- sigma21_k*Omega1rho1_k

      if (model=="S")
      {
         invAvec <- 1/(sigma22_k+vardir)                      
         invA  <- diag(invAvec)                                  
         first <- 1
         for (i in 1:D)
         {
            last <- first+Tmen1
            firstlast <- first:last
            invAZ1[firstlast,i]<-invAvec[firstlast]
            first <- first + T         
         }
      }
      else
      {
         rho2_k        <- thetak["rho2",1]                                      
         Unomenrho22_k <- 1-(rho2_k^2)                             
         Omega2drho2_k <- matrix(0,nrow=T,ncol=T)                       
         Omega2drho2_k[lower.tri(Omega2drho2_k)] <- rho2_k^seqTmen1_1
         Omega2drho2_k        <- Omega2drho2_k+t(Omega2drho2_k)    
         diag(Omega2drho2_k)  <- 1                                 
         Omega2drho2_k        <- (1/Unomenrho22_k)*Omega2drho2_k   
         sigma22Omega2drho2_k <- sigma22_k*Omega2drho2_k           
 
                                      #inverse in blocks
         first <- 1
         for (i in 1:D)
         {
            last  <- first+Tmen1
            firstlast <- first:last
            Ved   <- Ve[first:last,first:last]
            Ad    <- sigma22Omega2drho2_k + Ved
            invAd <- solve(Ad)              
            invA[first:last,first:last]<-invAd
            first <- first + T         
         }
         invAZ1 <- invA%*%Z1              
      }
      invVu1 <- solve(Vu1)

      invV <- invA - invAZ1%*%solve(invVu1+tZ1%*%invAZ1)%*%t(invAZ1)
      tXinvV     <- tX %*% invV
      inv_tXinVX <- solve(tXinvV %*% X) 
      P          <- invV - t(tXinvV) %*% inv_tXinVX %*% tXinvV

                                                # calculate S and F
      derivrho1_k <- p1derivrho1 + 2*rho1_k*tproxmatproxmat
                                                # calculate Va
      sigmaOmegaderivrho1Omega <- (-sigma21_k)*(Omega1rho1_k %*% 
                                       derivrho1_k %*% Omega1rho1_k)

      Va[[1]] <- Z1 %*% Omega1rho1_k %*% tZ1 
      Va[[2]] <- Z1 %*% sigmaOmegaderivrho1Omega %*% tZ1

      if (model=="S")
      {
         tZ1P    <- tZ1%*%P               
         tZ1PZ1  <- tZ1P%*%Z1           
         auxV1   <- tZ1PZ1%*%Omega1rho1_k       
         trPV[1] <- sum(diag(auxV1))   
         auxV2   <- tZ1PZ1%*%sigmaOmegaderivrho1Omega  
         trPV[2] <- sum(diag(auxV2))
         trPV[3] <- sum(diag(P))    
         Py      <- P %*% y                                    
         tyP     <- t(Py)                         # P is symmetric
         trPVPV[1,1] <- sum(auxV1*t(auxV1))
         trPVPV[2,2] <- sum(auxV2*t(auxV2))
         trPVPV[3,3] <- sum(P*t(P))        
         trPVPV[1,2] <- sum(auxV1*t(auxV2))
         tZ1PPZ1     <- tZ1P%*%t(tZ1P)         
         trPVPV[1,3] <- sum(tZ1PPZ1*t(Omega1rho1_k))
         trPVPV[2,3] <- sum(tZ1PPZ1*t(sigmaOmegaderivrho1Omega))
      }
      else
      {
         Va[[3]] <- diagonalizematrix(Omega2drho2_k,ntimes=D)

         derivOmega2drho2_k <- matrix(0,nrow=T,ncol=T)
         derivOmega2drho2_k[lower.tri(derivOmega2drho2_k)]<-seqTmen1_1*
                                               rho2_k^(seqTmen1_1-1) 
         derivOmega2drho2_k <- derivOmega2drho2_k +
                               t(derivOmega2drho2_k)
         derivOmega2drho2_k <- (1/Unomenrho22_k)*derivOmega2drho2_k +
                               (2*rho2_k/Unomenrho22_k)*Omega2drho2_k
         sigma22derivOmega2drho2_k <- sigma22_k * derivOmega2drho2_k

         Va[[4]] <- diagonalizematrix(sigma22derivOmega2drho2_k,ntimes=D)
         for (i in 1:nparam)
         {
            PV[[i]] <- P %*% Va[[i]]
            trPV[i] <- sum(diag(PV[[i]])) 
         }

         for (j in 1:nparam)
         {
            tPVj <- t(PV[[j]])
            for (i in 1:j)
               trPVPV[i,j] <- sum(PV[[i]]*tPVj)
         }
         Py  <- P %*% y  
         tyP <- t(Py)    
      }

      for (a in 1:nparam)
      {
         S[a] <- (-0.5)*trPV[a] + 0.5*(tyP %*% Va[[a]] %*% Py) 
         for (b in a:nparam)
            F[a,b] <- 0.5*trPVPV[a,b]
      }     

      for (a in 2:nparam)     # symmetric
         for (b in 1:(a-1))
            F[a,b] <- F[b,a]
                        
#      Finv <- try(solve(F), silent=TRUE)
      Finv <- try(ginv(F))
      if(class(Finv)=="try-error")
      {
         result$fit$convergence <- FALSE
         return (result)
      }
                     
      thetakmas1 <- thetak + Finv %*% S 

      # Restricting the spatial correlation to (-0.999,0.999)     
      if (thetakmas1["rho1",1]<= -1 )
         thetakmas1["rho1",1] <- -0.999
      else if (thetakmas1["rho1",1]>= 1 )
         thetakmas1["rho1",1] <- 0.999

      if (model=="ST")
      {  
         if (thetakmas1["rho2",1]<= -1 )
            thetakmas1["rho2",1] <- -0.999
         else if (thetakmas1["rho2",1]>= 1 )
            thetakmas1["rho2",1] <- 0.999
      }

      # Test values!=0 to avoid division errors
      if (any(thetak==0))
         for (i in 1:nparam)
            if (thetak[i]==0)
               thetak[i]<- 0.0001
      diff <- max( abs((thetak - thetakmas1)/thetak ))

   } #while (diff>PRECISION & k<MAXITER)

   thetakmas1["sigma21",1] <- max(thetakmas1["sigma21",1],0)
   thetakmas1["sigma22",1] <- max(thetakmas1["sigma22",1],0)
   sigma21_k <- thetakmas1["sigma21",1]
   sigma22_k <- thetakmas1["sigma22",1]

   if (thetakmas1["rho1",1]==-0.999)
      thetakmas1["rho1",1] <- -1
   else if (thetakmas1["rho1",1]==0.999)
      thetakmas1["rho1",1] <- 1
   rho1_k    <- thetakmas1["rho1",1]

   if (model=="ST")
   {
      if (thetakmas1["rho2",1]==-0.999)
         thetakmas1["rho2",1] <- -1
      else if (thetakmas1["rho2",1]==0.999)
         thetakmas1["rho2",1] <- 1
      rho2_k  <- thetakmas1["rho2",1]
   }


   # Indicator of convergence
   result$fit$iterations <- k
   if (k>=MAXITER && diff>=PRECISION)
   {
      result$fit$convergence <- FALSE
      return (result)
   }

   # Este if se supone que ya no hace falta 
   if (sigma21_k<0 || rho1_k < (-1) || rho1_k>1 || sigma22_k<0 || 
       (model=="ST" && (rho2_k<(-1) || rho2_k>1)) )  
   {
      result$fit$estvarcomp  <- data.frame(estimate=thetakmas1,std.error=0)
      print("Este mensaje no deberia salir... y ya no sale")
print(result$fit$estvarcomp)
      return(result)
   }
                                   # calculate estimates beta, u, mu
   Omega1rho1_k <- try(solve(crossprod(Id-rho1_k*proxmat))) 
   if(class(Omega1rho1_k)=="try-error")
   {
      result$fit$convergence <- FALSE
      result$fit$estvarcomp  <- data.frame(estimate=thetakmas1,std.error=0)
      return (result)
   }
   
   if (model=="S")
   {
      invAvec <- 1/(sigma22_k+vardir)
      invA <- diag(invAvec)            
      first <- 1
      for (i in 1:D)
      {
         last <- first+Tmen1
         firstlast <- first:last
         invAZ1[firstlast,i]<-invAvec[firstlast]
         first <- first + T         
      }
   }
   else
   {
#      rho2_k        <- thetakmas1[4] 
      Unomenrho22_k <- 1-(rho2_k^2)  
      Omega2drho2_k <- matrix(0,nrow=T,ncol=T)
      Omega2drho2_k[lower.tri(Omega2drho2_k)] <- rho2_k^seqTmen1_1
      Omega2drho2_k        <- Omega2drho2_k+t(Omega2drho2_k)      
      diag(Omega2drho2_k)  <- 1                              
      Omega2drho2_k        <- (1/Unomenrho22_k)*Omega2drho2_k
      sigma22Omega2drho2_k <- sigma22_k*Omega2drho2_k        
 
      first <- 1
      for (i in 1:D)
      {
         last <- first+Tmen1
         firstlast <- first:last
         Ad <- sigma22Omega2drho2_k + Ve[first:last,first:last]
         invAd <- solve(Ad)                                    
         invA[first:last,first:last]<-invAd                    
         first <- first + T         
      }
      invAZ1 <- invA%*%Z1
   }

   Vu1 <- sigma21_k*Omega1rho1_k
   if (sigma21_k!=0)
   {
      invVu1  <- try(solve(Vu1),silent=TRUE)    ####
      if(class(invVu1)=="try-error")
      {
         result$fit$convergence <- FALSE
         result$fit$estvarcomp  <- data.frame(estimate=thetakmas1,std.error=0)
      print ("se sale por la inversa de invVu1")
         return (result)
      }
      invV    <- invA - invAZ1%*%solve(invVu1+tZ1%*%invAZ1)%*%t(invAZ1)
   }
   else
      invV <- invA

   tXinvV  <- tX %*% invV
   Q       <- solve(tXinvV %*% X)                  
   betaest <- Q %*% (tXinvV %*% y)             
   ymenXbetaest  <- ( y - X %*% betaest )
   invVymenXBest <- invV %*% ymenXbetaest
   
   parte1  <- Vu1 %*% tZ1
   if (model=="S")
      parte2 <- diag(sigma22_k,nrow=M)
   else
      parte2 <- diagonalizematrix(sigma22Omega2drho2_k,ntimes=D)

   u1est   <- parte1 %*% invVymenXBest
   u2dtest <- parte2 %*% invVymenXBest
   u1dtest <- matrix(data=rep(u1est, each=T),nrow=M,ncol=1) 
   mudtest <- X%*%betaest + u1dtest + u2dtest                   
 
   V       <- solve(invV)
   loglike <- (-0.5) * ( M*log(2*pi) +
              determinant(V,logarithm=TRUE)$modulus +
              t(ymenXbetaest)%*%invV%*%ymenXbetaest )
   AIC     <- (-2)*loglike + 2*(length(betaest)+nparam)
   BIC     <- (-2)*loglike + log(M)*(length(betaest)+nparam)

                        # calculate pvalues
   std.errortheta <- sqrt(diag(Finv))
   std.errorbeta  <- sqrt(diag(Q))
   tvalue         <- betaest/std.errorbeta
   pvalue         <- 2*pnorm(abs(tvalue),lower.tail=FALSE)
   
   result$fit$estvarcomp <- data.frame(estimate=thetakmas1,
                                     std.error=std.errortheta)
   result$fit$estcoef    <- data.frame(beta=betaest,
                                     std.error=std.errorbeta,
                                     tvalue=tvalue,
                                     pvalue=pvalue)
   result$fit$goodness   <- c(loglike=loglike, AIC=AIC, BIC=BIC)
   result$eblup          <- mudtest
   return (result)
}

