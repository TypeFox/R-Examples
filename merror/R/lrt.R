"lrt" <-
function(x,M=40) 
{ 
# Likelihood-ratio test for beta.bars for NonConstant Bias model using MLE
# Jaech, pp. 204-205
# sigma2.cb and sigma.mu2.cb and sigma2.ncb and sigma.mu2.ncb should be the MLE

  N <- dim(x)[2]
  n <- dim(x)[1]
  sig <- matrix(NA,N,N)
  beta.bars <- beta.bar(x)
  
  ncb.model <- ncb.od(x,M=M)
  cb.model  <- ncb.od(x,beta=rep(1,N),M=M)

  sigma2.ncb <- ncb.model$sigma2
  sigma.mu2.ncb <- ncb.model$sigma.mu2  
  
  sigma2.cb <- cb.model$sigma2
  sigma.mu2.cb <- cb.model$sigma.mu2  
 
# cat("\nbeta.bars=",beta.bars)
# cat("\nsigma.mu2.ncb=",sigma.mu2.ncb)
# cat("\nsigma.mu2.cb=",sigma.mu2.cb)
# cat("\nsigma2.ncb=",sigma2.ncb)
# cat("\nsigma2.cb=",sigma2.cb)

  # Likelihood L.mle for ncb mle
  
  a <- beta.bars*sqrt(sigma.mu2.ncb)
  
# cat("\na=",a)

  V <- prod(sigma2.ncb)
  Q <- V*(sum(a^2/sigma2.ncb)+1)
  
# cat("\nV=",V)
# cat("\nQ=",Q)

  for(i in 1:N)
  {
    for(j in 1:N)
    {
      if(j==i) sig[i,j] <- V*(sum(a[-i]^2/sigma2.ncb[-i])+1)/(Q*sigma2.ncb[i])
      else sig[i,j] <- -V*a[i]*a[j]/(Q*sigma2.ncb[i]*sigma2.ncb[j])
    }
  }
  
# cat("\nsig\n")
# print(sig)
# cat("\n")
  
  L.mle <- -0.5*n*log(Q) - 0.5*(n-1)*sum(sig*var(x))
  
# cat("\n-0.5*n*log(Q)=",-0.5*n*log(Q))
# cat("\n-0.5*(n-1)*sum(sig*var(x))=",-0.5*(n-1)*sum(sig*var(x)))
  
# cat("\nL.mle=",L.mle)

  # Likelihood L.hyp for constant bias

  a <- rep(1,N)*sqrt(sigma.mu2.cb)

# cat("\na=",a)

  V <- prod(sigma2.cb)
  Q <- V*(sum(a^2/sigma2.cb)+1)
  
# cat("\nV=",V)
# cat("\nQ=",Q)

  for(i in 1:N)
  {
    for(j in 1:N)
    {
      if(j==i) sig[i,j] <- V*(sum(a[-i]^2/sigma2.cb[-i])+1)/(Q*sigma2.cb[i])
      else sig[i,j] <- -V*a[i]*a[j]/(Q*sigma2.cb[i]*sigma2.cb[j])
    }
  }
  
# cat("\nsig\n")
# print(sig)
# cat("\n")
  
  L.hyp <- -0.5*n*log(Q) - 0.5*(n-1)*sum(sig*var(x))
  
# cat("\n-0.5*n*log(Q)=",-0.5*n*log(Q))
# cat("\n-0.5*(n-1)*sum(sig*var(x))=",-0.5*(n-1)*sum(sig*var(x)))
  
# cat("\nL.hyp=",L.hyp,"\n\n")
  
  lambda <- -2*(L.hyp - L.mle)
  
# cat("\nlambda = ",lambda,"\n")
  
  # df for Likelihood-Ratio Test = no. of parameters lost
  # For this test there are 2*N unrestricted parameters 
  #   (N sigma2s and N betas) and
  #   there are only N + 1 parameters in the restricted model
  #   because all the betas are reduced to just one. Thus
  #   N - 1 parameters are lost (or 2*N - (N + 1) = 2*N - N - 1 = N - 1)
  # So the degrees of freedom = N - 1
  
  df <- N - 1
  
# cat("\ndf = ",df,"\n")
  
  list(N.methods=N,n.items=n,beta.bars=beta.bars,lambda=lambda, df=df,
  p.value=1-pchisq(lambda,df)) 
}
