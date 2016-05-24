
 
pred.dist.simul = function(hyperest, tpt, include.obs = T, N =1)
{
#
# This function simulates N- replicates from the predictive distribution for
# a given time point (tpt) from 1 to n (length of the data).
# 
# Input
#   hyperest: Output from the "staircase.hyper.est" functions, containing 
#             estimates of all hyperparameters
#   tpt:      A specific time point - from 1 to n corresponding to the 
#             number of time points from the data set 
#   include.obs:  If True, the observed data, for time "tpt", are also returned
#   N:        Number of replicates

# Output:
#   A matrix with N rows; the number of columns depends on whether the observed data are returned 
#   The columns are organized consistent with the observed data 
#         (ie. uxp ungauged blocks, g1xp, g2xp , ...)

# Note: This function could be slow if there are missing data at gauged sites 
#       correspondind to the selected time point. That is, it is fastest at time points
#       corresponding to Block 1 and slower with higher blocks.


y = as.matrix(hyperest$data)
b0 = hyperest$Beta0
block = hyperest$block
delta = hyperest$Delta
lambda = hyperest$Lambda
omega = hyperest$Omega
ZZ = hyperest$covariate

n = dim(y)[1]
p = dim(hyperest$Omega)[1]
m = unique(apply(is.na(y),2,sum))
gm = c(0, cumsum(hyperest$block*p))
K = length(m)

# Identify the proper block
blk = K
if (tpt > m[1]) blk = 0
if (blk > 0) while ((tpt > m[blk])&(blk > 0))  blk = blk -1
  
# When tpt > m[1] - ie. no missing data at any steps.                                                                  
 if (blk == 0) y.g = matrix(rep(y[tpt,],N),byrow=T, nrow = N)

 if (blk > 0) {
# Generate data for the gauged locations
   y.g = NULL
   for (i in 1:N) {
    # Start with last block
    if (m[K] >0) {
      # generate independent realizations from t-distribution with d0 degree of freedom
       d0 = delta[[K]]-block[K]*p + n-m[K]+1
       resid = matrix(rt(m[K]*block[K]*p, d0), nrow=m[K])
         # getting predictive mean and covariance
           y1 = y[(m[K]+1):n,(gm[K]+1):(gm[K+1])]
           mu = matrix(0,nrow=n,ncol= gm[K+1] - gm[K])
           if (!is.null(ZZ)) mu = ZZ %*% b0[,(gm[K]+1):(gm[K+1])]
            mu1 = (mu)[1:(m[K]),]
            mu2 = (mu)[(m[K]+1):n,]
         
           if (!is.null(ZZ)) A = diag(n) + ZZ %*% hyperest$Finv %*% t(ZZ)
               else A = diag(n)
               A11 = A[1:(m[K]),1:(m[K])]
               A12 = A[1:(m[K]),(m[K]+1):n]
               A22 = A[(m[K]+1):n,(m[K]+1):n]
            phi.K = (A11 - A12 %*% solve(A22) %*% t(A12))* (delta[[K]] - block[K]*p +1)/d0
            Psi.K = (kronecker(lambda[[K]],omega)+ t(y1 - mu2) %*% solve(A22) %*% (y1 - mu2))/
                         (delta[[K]]-block[K]*p +1)
            mu.K = mu1 + A12 %*% solve(A22) %*% (y1 - mu2)
 
     Tb = chol(Psi.K)
     Ta = chol(phi.K)
     y[1:m[K],(gm[K]+1):(gm[K+1])] = mu.K + t(Ta) %*% resid %*% Tb
     } 
    # The remaining blocks
    for (j in (K-1):1) {
          d0 = delta[[j]]-block[j]*p + n-m[j]+1
          resid = matrix(rt(m[j]*block[j]*p, d0), nrow=m[j])
         # getting predictive mean and covariance
           y1 = y[(m[j]+1):n,(gm[j]+1):(gm[j+1])]
           mu = matrix(0,nrow=n,ncol= gm[j+1] - gm[j])
           if (!is.null(ZZ)) mu = ZZ %*% b0[,(gm[j]+1):(gm[j+1])]
            mu1 = (mu)[1:(m[j]),]
            mu2 = (mu)[(m[j]+1):n,]

           mu.1 = matrix(0,nrow=n,ncol= gm[K+1] - gm[j+1])
           if (!is.null(ZZ)) mu.1 = ZZ %*% b0[,(gm[j+1]+1):(gm[K+1])]
           eps.t = y[,(gm[j+1]+1):(gm[K+1])] - mu.1
         
           if (!is.null(ZZ)) A = diag(n) + ZZ %*% hyperest$Finv %*% t(ZZ) +
                                 eps.t %*% solve(hyperest$Hinv[[j]]) %*% t(eps.t)
               else A = diag(n)
               A11 = A[1:(m[j]),1:(m[j])]
               A12 = A[1:(m[j]),(m[j]+1):n]
               A22 = A[(m[j]+1):n,(m[j]+1):n]
            phi.j = (A11 - A12 %*% solve(A22) %*% t(A12))* (delta[[j]] - block[j]*p +1)/d0
            Psi.j = (kronecker(lambda[[j]],omega)+ t(y1 - mu2) %*% solve(A22) %*% (y1 - mu2))/
                         (delta[[j]]-block[j]*p +1)
            mu.j = mu1 + A12 %*% solve(A22) %*% (y1 - mu2)
 
     Tb = chol(Psi.j)
     Ta = chol(phi.j)
     y[1:m[j],(gm[j]+1):(gm[j+1])] = mu.j + t(Ta) %*% resid %*% Tb
     } 
    
y.g = rbind(y.g,y[tpt,])

}
}

# Generate observations at ungauged sites conditional on y.g

 delta.0 = hyperest$Delta.0
 u = dim(hyperest$Lambda.0)[1]
 d0.0 = delta.0 - u*p+1
 resid = matrix(rt(N*u*p, d0.0), nrow=N)
 b0.0 = matrix(rep(b0[,1:p],u),nrow=dim(b0[,1:p])[1])
 Zt = ZZ[tpt,]
 ZZ.0 = matrix(rep(Zt, N), byrow=T, ncol=length(Zt))
 mu.0 = ZZ.0 %*% b0.0 + (y.g - ZZ.0 %*% b0) %*% kronecker(hyperest$Xi0.0, diag(p))
 
 Psi.0 = (kronecker(hyperest$Lambda.0,omega))/(d0.0)
 phi.0 =  diag(N) + ZZ.0 %*%  hyperest$Finv %*% t(ZZ.0) + 
               (y.g - ZZ.0 %*% b0) %*% hyperest$H.0 %*% t(y.g - ZZ.0 %*% b0)

 Tb = chol(Psi.0)
 if (dim(phi.0)[1] > 1) Ta = chol(diag(diag(phi.0)))
 if (dim(phi.0)[1] == 1) Ta = sqrt(phi.0)

 y.u = mu.0 + t(Ta) %*% resid %*% Tb

obj = cbind(y.u,y.g)
if (!include.obs) obj = cbind(y.u,y.g[,1:(gm[j+1])])
return(obj)
}



