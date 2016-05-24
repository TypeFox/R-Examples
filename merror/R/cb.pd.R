"cb.pd" <-
function(x, conf.level=0.95, M=40)
{
# Maximum Likelihood precision estimates for Constant Bias Model using Paired Data
# ME (method of moments estimator) and MLE are the same for i=3 instruments
#  except for a factor of (n-1)/n: MLE = (n-1)/n * ME
#
# Using paired differences forces Constant Bias model
#
# This function first computes the Grubbs ME and then uses it as the starting
#   point for iteratively computing the MLEs (required for N > 3)
# 
#
# Jaech 1985, Chapters 3, 4, and 5 p. 144 in particular, eq. 5.3.1
#
# Written for any no. of methods N
#
# x[i,k] = alpha[i] + beta[i]*mu[k] + epsilon[i,k]
#   with beta[1] = beta[2] = ... = beta[N]
#   N = no. of methods
#   n = no. of items
# 

  n <- dim(x)[1]
  N <- dim(x)[2]

#  Compute alphas for unconstrained betas (NCB) and betas=1 (CB) cases
#    for completeness 

   alpha.cb  <- colMeans(x) - rep(1,N)*mean(colMeans(x))
   alpha.ncb <- colMeans(x) - beta.bar(x)*mean(colMeans(x))

# compute paired differences 

  y <- array(NA,c(N,N,n)) # not all elements are used
  
  for(i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      y[i,j,] <- x[,i] - x[,j]
    }
  }
  
# print(y[1,2,])

# compute variances for each paired difference

  V <- matrix(NA,N,N) # not all elements are used
  
  for(i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      V[i,j] <- var(y[i,j,])
      V[j,i] <- var(y[i,j,])
    }
  }
  
#  print(V)

# Compute Grubbs (Squared) Precions estimates - stored in r

  S <- vector("numeric",N)
  
  for(i in 1:N)
  {
    S[i] <- 0
    for(j in 1:N)
    {
      if(j!=i) S[i] <- S[i] + V[i,j]
    }
  }
  
# print(S)
  
  VT <- sum(S)/2  

  r <- vector("numeric",N)
  
  for(i in 1:N)
  {
    r[i] <- ((N-1)*S[i] - VT)/((N-1)*(N-2))
  }

# cat("\nGrubbs initial r =",r)

  initial.r <- r

# Use Grubbs MEs for squared precisions to find MLEs interatively

# compute MLE

  for(i in 1:M) r <- mle(V, r, n)
# cat("\nr =",r)

# compute squared se's for MLE
  sigma2.se2 <- mle.se2(V,r,n)
# cat("\nsigma2.se2 =",sigma2.se2)

# compute degrees of freedom

  df <- 2*r^2/sigma2.se2
  
# cat("\ndf =",df,"\n")

# compute confidence intervals

  sig.level <- 1 - conf.level

  lb <- df*r/qchisq(1-sig.level/2,df)
  ub <- df*r/qchisq(sig.level/2,df)

  sigma.table <- data.frame(n=rep(n,N),sigma=sqrt(r),sigma.se=sigma2.se2^(1/4),
    alpha.cb=alpha.cb,alpha.ncb=alpha.ncb,
    beta=rep(1,N),df=df,chisq.low=qchisq(sig.level/2,df),
    chisq.upper=qchisq(1-sig.level/2,df),lb=sqrt(lb),ub=sqrt(ub))
      
  dimnames(sigma.table)[[1]] <- names(x)

# return all the table and all the pieces

 
list(
    conf.level=conf.level,
    sigma.table=sigma.table,
    n.items=n,
    N.methods=N,
    Grubbs.initial.sigma2=initial.r,
    sigma2=r,
    sigma2.se2=sigma2.se2,
    alpha.cb=alpha.cb,
    alpha.ncb=alpha.ncb,
    beta=rep(1,N),df=df,
    chisq.low=qchisq(sig.level/2,df),
    chisq.upper=qchisq(1-sig.level/2,df),
    lb=lb,
    ub=ub
    )
  
}
