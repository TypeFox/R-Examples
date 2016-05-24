"precision.mle.ncb.od" <-
function(x, M=20, beta.bars=beta.bar(x),jaech.errors=FALSE) 
{
# Iterative Approximation to MLE precision estimates for NonConstant Bias model
#   using Original Data
# Jaech, p. 185-186
# Use Grubbs NonConstant Bias using Original Data for initial values for
#   precision (sigma)
# Default for beta.bars is NonConstant Bians (betas differ from 1.0) by
#   using Grubbs least squares estimate via function "beta.bar"
# To compute under assumption of constant bias, use a vector of N ones.

  N <- dim(x)[2]  # no. of instruments or methods
  n <- dim(x)[1]  # no. of items
  
# cat("\nPrecision Estimates Using MLE\nAssuming A Model With NonConstant Bias Using Original Data\n")
# cat("Jaech, Chapter 6, p. 185\n")
# cat("Note errors in Jaech when estimating sigma.mu^2 - See function process.var.mle.jaech.err\n")
  
  if(jaech.errors==TRUE)
    cat("\n***Using Same Errors in Jaech's Fortran Program (p. 288) for sigma.mu^2 for compatibility***\n")
  
# cat("\nNo. of methods N=",N,"\nNo. of items n=",n,"\nNo. of iterations for MLE M=",M,"\n")

  s <- var(x)
# cat("\nvariance-covariance matrix s\n")
# print(s)
# cat("\n")

  sigma2 <- round(precision.grubbs.cb.pd(x),4) # initial precision estimates 1 to N
# sigma2[sigma2<0] <- abs(sigma2[sigma2<0]) # take absolute value (remove negative values)
# cat("\nInitial Squared Precision Estiamtes (Grubbs)\n")
# print(sigma2)
# cat("\n")
  
# cat("\nDefault Least Squares Estiamtes of Betas (Grubbs) -or- Assumed Values\n")
# print(beta.bars)
# cat("\n")
    
  if(jaech.errors==TRUE) sigma.mu2 <- process.var.mle.jaech.err(sigma2,s,beta.bars,N,n)
  else sigma.mu2 <- process.var.mle(sigma2,s,beta.bars,N,n)
  
# if(sigma.mu2 < 1) sigma.mu2 <- 0 # set sigma.mu2 = 0 if negative
  
# cat("\nInitial Squared Precision Estimates (Grubbs)\n")
# cat(0,sigma2,"sigma.mu^2=",sigma.mu2,"\n")
# cat("\n")
  
  for(m in 1:M)
  {
    for(i in 1:N)
    {
      sigma2[i] <- sigma_mle(i,s,sigma2,sigma.mu2,beta.bars,N,n)
#      cat("\n","i=",i,"sigma2",sigma2)
    }
#   cat("\n",m,"Sq.Prec.=",sigma2,"sigma.mu^2=",sigma.mu2)
        
    if(jaech.errors==TRUE) sigma.mu2 <- process.var.mle.jaech.err(sigma2,s,beta.bars,N,n)
    else sigma.mu2 <- process.var.mle(sigma2,s,beta.bars,N,n)
  }
# cat("\n\n")    
# cat("\nPrecision (Final)=",sqrt(sigma2),"sigma.mu (Final)=",sqrt(sigma.mu2),"\n\n")

  list(sigma2=sigma2,sigma.mu2=sigma.mu2)

}
