pbmseSTFH <-
function(formula,D,T,vardir,proxmat,B=100,model="ST",MAXITER=100,PRECISION=0.0001,data)
{
   result <- list(est=NA, mse=NA)

   if (!is.matrix(proxmat))
      proxmat <- as.matrix(proxmat)

   namevar <- deparse(substitute(vardir))
   if (!missing(data))
   {
      formuladata <- model.frame(formula,na.action = na.omit,data)
      Xdt           <- model.matrix(formula,data)        
      vardir    <- data[,namevar]
   } else
   {
      formuladata <- model.frame(formula,na.action = na.omit)
      Xdt           <- model.matrix(formula)        
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

   nformula  <- nrow(Xdt)   
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


   result$est <- eblupSTFH(y~Xdt-1,D,T,vardir,proxmat,model,MAXITER,PRECISION)
   if (result$est$fit$convergence==FALSE) 
   {
      warning("The fitting method does not converge.\n") 
      return (result);
   }
   beta       <- result$est$fit$estcoef[,"beta"]
   estvarcomp <- result$est$fit$estvarcomp[,"estimate"]

   M <- D*T
   sigma21 <- estvarcomp[1]
   rho1    <- estvarcomp[2]
   sigma22 <- estvarcomp[3]

   msedt <- 0
   Id    <- diag(1,nrow=D, ncol=D)
   Omega1rho1 <- try(solve(crossprod(Id-rho1*proxmat))) 
   if (class(Omega1rho1)=="try-error")
   {
      result$mse <- 0
      return (result)
   }
   sigma21Omega1rho1 <- sigma21*Omega1rho1
 
   if (model=="ST")
   {
      rho2   <- estvarcomp[4]
      Unomenrho22_05 <- (1-rho2^2)^(-0.5)
      u2dt_b <- matrix(0, nrow=M, ncol=1)
   }

   cat("\nBootstrap procedure with B =",B,"iterations starts.\n")

   b<-1
   while (b<=B)
   {
      u1_b   <- matrix(data=mvrnorm(n=1, mu=rep(0,D),
                       Sigma=sigma21Omega1rho1), nrow=D, ncol=1)
      u1dt_b <- matrix(data=rep(u1_b, each=T),nrow=M,ncol=1)   

      if (model=="S")       
         u2dt_b <- matrix(data=rnorm(M, mean=0, sd=sqrt(sigma22)), 
                          nrow=M, ncol=1) 
      else
      {
         epsilondt_b <- matrix(data=rnorm(M, mean=0, sd=sqrt(sigma22)),
                               nrow=M, ncol=1)
         i <- 1
         for (d in 1:D)
         {
            u2dt_b[i] <- Unomenrho22_05*epsilondt_b[i]
            for (t in 2:T)
            {
               i <- i+1
               u2dt_b[i] <- rho2*u2dt_b[i-1]+epsilondt_b[i]
            }
            i<-i+1
         }
      }

      edt_b  <- matrix(data=rnorm(M, mean=0, sd=sqrt(vardir)),
                       nrow=M,ncol=1) 
      ydt_b  <- Xdt%*%beta + u1dt_b + u2dt_b + edt_b
      mudt_b <- ydt_b - edt_b
                                               ### fitting the model  
      resultb <- eblupSTFH(ydt_b~Xdt-1,D,T,vardir,proxmat,model,MAXITER,PRECISION) 
 
      if (resultb$fit$convergence==FALSE)
         next

         # Test estimated components of variance are valids
      if (resultb$fit$estvarcomp["sigma21","estimate"]<0 | 
          resultb$fit$estvarcomp["rho1","estimate"]<(-1) | 
          resultb$fit$estvarcomp["rho1","estimate"]>1    | 
          resultb$fit$estvarcomp["sigma22","estimate"]<0 | 
          (model=="ST" & 
           (resultb$fit$estvarcomp["rho2","estimate"]<(-1) | 
            resultb$fit$estvarcomp["rho2","estimate"]>1)
          ))    
      {
         print (resultb$fit$estvarcomp)
         cat("pbmseSTFH: Este  mensaje no deberia salir: fuera del rango\n")
         next
      }
      cat("b =",b,"\n")     
                                                ### calculate msedt
      difference <- resultb$eblup - mudt_b
      msedt      <- msedt + difference^2
      b <- b+1    
   } #while (b<=B)

   result$mse <- msedt/B
   return (result)
}

