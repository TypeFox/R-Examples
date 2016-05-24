mseFH <-
function(formula,vardir,method="REML",MAXITER=100,PRECISION=0.0001,B=0,data)
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


   result$est <- eblupFH(y~X-1,vardir,method,MAXITER,PRECISION,B)
   if (result$est$fit$convergence==FALSE) 
   {
      warning("The fitting method does not converge.\n") 
      return (result);
   }

   A  <- result$est$fit$refvar

   m<-dim(X)[1] # Sample size or number of areas
   p<-dim(X)[2] # Num. of auxiliary variables (including intercept)

  # Initialize vectors containing the values of g1-g3 for each area along with the final mse 
   g1d<-rep(0,m)
   g2d<-rep(0,m)
   g3d<-rep(0,m)
   mse2d<-rep(0,m)
  
  # Elements of the inverse covariance matrix in a vector
   Vi<-1/(A+vardir)
  # Auxiliary calculations
   Bd<-vardir/(A+vardir)
   SumAD2<-sum(Vi^2)
   XtVi<-t(Vi*X)
   Q<-solve(XtVi%*%X)
   

  # Calculation of g1-g3 and final MSE when fitting method is REML
   if (method=="REML")
   {
      VarA<-2/SumAD2      # Asymptotic variance of REML/ML estimator of variance A
      for (d in 1:m){
         g1d[d]<-vardir[d]*(1-Bd[d])
         xd<-matrix(X[d,],nrow=1,ncol=p)
         g2d[d]<-(Bd[d]^2)*xd%*%Q%*%t(xd)
         g3d[d]<-(Bd[d]^2)*VarA/(A+vardir[d])
         mse2d[d]<-g1d[d]+g2d[d]+2*g3d[d]
      }
    
   } else if (method=="ML")
   {
      VarA<-2/SumAD2      # Asymptotic variance of REML/ML estimator of variance A
      b   <-(-1)*sum( diag( Q%*%(t((Vi^2)*X)%*%X) ) )/SumAD2
      for (d in 1:m){
         g1d[d]<-vardir[d]*(1-Bd[d])
         xd<-matrix(X[d,],nrow=1,ncol=p)
         g2d[d]<-(Bd[d]^2)*xd%*%Q%*%t(xd)
         g3d[d]<-(Bd[d]^2)*VarA/(A+vardir[d])
         mse2d[d]<-g1d[d]+g2d[d]+2*g3d[d]-b*(Bd[d]^2)
      }

   } else # Calculation of g1-g3 and final MSE when fitting method is FH
   {  
      SumAD<-sum(Vi)
      VarA<-2*m/(SumAD^2)               # Asymptotic variance of FH estimator of variance A
      b<-2*(m*SumAD2-SumAD^2)/(SumAD^3) # Asymptotic bias of FH estimator of A
      for (d in 1:m){  
         g1d[d]<-vardir[d]*(1-Bd[d])
         xd<-matrix(X[d,],nrow=1,ncol=p)                                              
         g2d[d]<-(Bd[d]^2)*xd%*%Q%*%t(xd)
         g3d[d]<-(Bd[d]^2)*VarA/(A+vardir[d])
         mse2d[d]<-g1d[d]+g2d[d]+2*g3d[d]-b*(Bd[d]^2)
       }
   } 

   result$mse <- mse2d   
   return (result)
}

