stoch.growth.rate<-function(matrices, prob=NULL, maxt=50000, verbose=TRUE)
{
   if(is.list(matrices)){matrices<-matrix(unlist(matrices), ncol=length(matrices))}
   s <- sqrt(dim(matrices)[1])   ## number of stage classes
   n  <-dim(matrices)[2]         ## number of matrixes
   #default equal probabilities
   if(is.null(prob)){prob<-rep(1/n, n)}
   Abar<-numeric(s^2)
   Exy <-numeric(s^4)
   ## for each matrix, add values to Abar and Exy weighted by probabilities in prob
   for (i in 1:n)
   {
      A<-matrices[,i]
      Exy <-Exy +prob[i] * kronecker(A,A)
      Abar<-Abar+prob[i] * A
   }
   ## Covariance matrix
   C<- (Exy - kronecker(Abar, Abar))*n/(n-1)
   C<-matrix(C, nrow=s^2)
   Abar<-matrix(Abar, nrow=s)
   ## code from eigen.analysis for lambda and Sensitivity matrix
   ev <- eigen(Abar)
    lmax <- which(Re(ev$values) == max(Re(ev$values)))
    lambda <- Re(ev$values[lmax])
    W <- ev$vectors
    w <- abs(Re(W[, lmax]))
    V <- Conj(solve(W))
    v <- abs(Re(V[lmax, ]))
    S <- v %o% w
   ##  Simulation
   r<-numeric(maxt)
   n0<-w
   for (t in 1:maxt)
   {
      if(verbose)
      {
         if(t==1 || t %% 10000 == 0){print(paste("Calculating stochastic growth at time", t), quote=FALSE)}
      }
      col<-sample(1:n, 1, prob=prob)
      A<-matrix(matrices[,col], nrow=s)
      n0<-A %*% n0
      N<-sum(n0)
      r[t]<-log(N)
      n0<-n0/N
   }   
   loglsim<-mean(r)
   dse<-1.96*sqrt(var(r)/maxt)
   CI<-c(loglsim-dse, loglsim+dse)
   ## Tuljapurkar approximation
   Svec<-matrix(S, ncol=1)
   tau2<-t(Svec) %*% C %*% Svec   
   loglams<-log(lambda) - tau2/(2*lambda^2)
   ## output...
   stoch<-list(
      approx=as.numeric(loglams),
      sim=loglsim,
      sim.CI=CI
               )
   stoch
}
