"ncb.od" <-
function(x, beta=beta.bar(x),M=40,conf.level=0.95)
{
# MLE Confidence Intervals for Precision
# Jaech 1985, pp. 201-202
#
# squared precisions stored in sigma2
# n is the no. of items
# N is the no. of methods
# beta is a vector of slope biases

   n <- dim(x)[1]
   N <- dim(x)[2]
   s <- var(x)

   H <- matrix(0,N+1,N+1)

#  Compute alphas for unconstrained betas (NCB) and betas=1 (CB) cases
#    for completeness

   alpha.cb  <- colMeans(x) - rep(1,N)*mean(colMeans(x))
   alpha.ncb <- colMeans(x) - beta.bar(x)*mean(colMeans(x))


   mles <- precision.mle.ncb.od(x, beta.bars=beta,M=M)

   sigma2 <- mles[[1]]
#  cat("\nsigma2 =",sigma2)
   sigma.mu2 <- mles[[2]]
#  cat("\nsigma.mu2 =",sigma.mu2)

   a <- beta*sqrt(sigma.mu2)
   c0 <- a^2
   d2 <- sum(beta^2/sigma2)
   c1 <- vector("numeric",N)
   c2 <- vector("numeric",N)
   c3 <- vector("numeric",N)
   c4 <- vector("numeric",N)
   b0 <- vector("numeric",N)
   b1 <- vector("numeric",N)
   b2 <- vector("numeric",N)


   for(i in 1:N)
   {
  #  cat("\ni=",i)
     c1[i] <- sum(a[-i]^2/sigma2[-i]) + 1
  #  cat("\nc1=",c1[i])
     c2[i] <- sum(diag(s[-i,-i])/sigma2[-i])
  #  cat("\nc2=",c2[i])
     c3[i] <- sum(((c1[i] - a[-i]^2/sigma2[-i])*diag(s[-i,-i]))/sigma2[-i])
  #  cat("\nc3=",c3[i])
     c4[i] <- sum(a[-i]*s[i,-i]/sigma2[-i])
  #  cat("\nc4=",c4[i])
     b0[i] = c0[i]
  #  cat("\nb0=",b0[i])
     b1[i] = c1[i]
  #  cat("\nb1=",b1)
     b2[i] = c1[i]*s[i,i] + c0[i]*c2[i] - 2*a[i]*c4[i]
  #  cat("\nb2=",b2[i])

     H[i,i] <- -0.5*n*b1[i]^2/(b0[i] + b1[i]*sigma2[i])^2
  #  cat("\nH[",i,",",i,"]=",H[i,i])
   }

   for(i in 1:(N-1))
   {
     for(j in (i+1):N)
     {
    #  cat("\nj=",j)
       g1 <- -beta[j]^2*sigma.mu2/sigma2[j]^2
    #  cat("\ng1=",g1)
       g2 <- sigma.mu2*(-beta[j]^2*s[i,i] - beta[i]^2*s[j,j] + 2*beta[i]*beta[j]*s[i,j])/sigma2[j]^2
    #  cat("\ng2=",g2)

       g3.b <- 0

       for(k in 1:N)
       {
         if(k!=i&k!=j) g3.b <- g3.b + (beta[k]^2*s[j,j] + beta[j]^2*s[k,k] - 2*beta[j]*beta[k]*s[j,k])/sigma2[k]
       }

       g3 <- -s[j,j]/sigma2[j]^2 - sigma.mu2/sigma2[j]^2*g3.b
    #  cat("\ng3=",g3,"\n\n")

       H[i,j] <- n*beta[j]^2*sigma.mu2*(b1[i]*sigma2[i] + 0.5*b0[i])/(sigma2[j]^2*(b0[i] + b1[i]*sigma2[i])^2) -
   0.5*(n-1)*(b0[i]*g3 - b1[i]*g2 - b2[i]*g1)/(b0[i] + b1[i]*sigma2[i])^2
       H[j,i] <- H[i,j]
    #  cat("\nH[",i,",",j,"]=",H[i,j])
     }
   }

   H[N+1,N+1] <- -0.5*n*d2^2/(d2*sigma.mu2 + 1)^2
#  cat("\nH[",N+1,",",N+1,"]=",H[N+1,N+1])

   for(i in 1:N)
   {
     H[i,N+1] <- (n*beta[i]^2*(d2*sigma.mu2 + 0.5) - (n - 1)*beta[i]*(beta[i]*s[i,i]/sigma2[i] + sum(beta[-i]*s[i,-i]/sigma2[-i])))/
       (sigma2[i]^2*(d2*sigma.mu2 + 1)^2)
     H[N+1,i] <- H[i,N+1]
  #  cat("\nH[",i,",",N+1,"]=",H[i,N+1])
   }
   H <- -H

#  cat("\nH matrix\n")
#  print(H)
#  cat("\nVariances for Squared Precision Estimates\n")

   se2.sigma2 <- diag(solve(H))
   names(se2.sigma2) <- c(dimnames(x)[[2]],"Process")
   names(sigma2) <- dimnames(x)[[2]]
   names(beta) <- dimnames(x)[[2]]
#  print(se2.sigma2)

#  cat("\n\n")

   # Jaech 1985, p. 71, eq. 3.4.2
   df <- 2*c(sigma2^2,sigma.mu2^2)/se2.sigma2
   names(df) <- c(dimnames(x)[[2]],"Process")
   sig.level <- 1 - conf.level
   lb <- df*c(sigma2,sigma.mu2)/qchisq(1 - sig.level/2,df)
   ub <- df*c(sigma2,sigma.mu2)/qchisq(sig.level/2,df)

   out <- data.frame(n=rep(n,N+1),
     sigma=c(sqrt(sigma2),sqrt(sigma.mu2)),se.sigma=se2.sigma2^(1/4),
     alpha.cb=c(alpha.cb,NA),alpha.ncb=c(alpha.ncb,NA),beta=c(beta,NA),
     df,chisq.l=qchisq(sig.level/2,df),
     chisq.u=qchisq(1 - sig.level/2,df),lb=sqrt(lb),ub=sqrt(ub),
     bias.adj.sigma=c(sqrt(sigma2)/beta,NA))

   dimnames(out)[[1]] <- c(dimnames(x)[[2]],"Process")

# Compute errors
#   for no bias model (alpha=0, beta=1)

    errors.nb <- x - 1*apply(x,1,mean)

#   for constant bias model (alpha, beta=1)

    alpha.cb.mat <- matrix(alpha.cb,n,N,byrow=TRUE)
    errors.cb <- x - (alpha.cb.mat + 1*apply(x,1,mean))

#   for nonconstant bias model (alpha, beta)

    alpha.ncb.mat <- matrix(alpha.ncb,n,N,byrow=TRUE)
    errors.ncb <- x - (alpha.ncb.mat + matrix(apply(x,1,mean),n,N)*matrix(beta,n,N,byrow=TRUE))

   list(conf.level=conf.level,sigma.table=out,n.items=n,N.methods=N,sigma2=sigma2,sigma.mu2=sigma.mu2,
     se2.sigma2=se2.sigma2,alpha.cb=alpha.cb,alpha.ncb=alpha.ncb,
     beta=beta,df=df,lb=sqrt(lb),ub=sqrt(ub),
     bias.adj.sigma=c(sqrt(sigma2)/beta,NA),
     H=H,
     errors.nb=errors.nb,errors.cb=errors.cb,errors.ncb=errors.ncb)

}
