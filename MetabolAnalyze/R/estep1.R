estep1 <-
function(Y, Tau, Pi, mu, W, Sig, g, p, reset)
{
  logTau<-Tau
  for(i in 1:g)
  {
     if(g==1)
     {
     	Tau<-matrix(rep(1,nrow(Y)))
      }else{
        Tau[,i]<-(dmvnorm(Y, mu[,i], W[,,i]%*%t(W[,,i])+ Sig*diag(p)))*Pi[i]
        logTau[,i]<-(dmvnorm(Y, mu[,i], W[,,i]%*%t(W[,,i])+ Sig*diag(p), log=TRUE)) + log(Pi[i])
      }
  }
  Tau[,i][Tau[,i]==Inf] <- 1
  Tau[,i][Tau[,i]==-Inf] <- 0
  Tau<-Tau/apply(Tau,1,sum)
  ## In case the posterior probability of membership (Tau) is NA due to some computational problems 
  ## then assign the observation into the group in which it has highest log of Tau.
  if(sum(is.na(Tau)) != 0) 
  {                                                                                       
  	reset<-TRUE
  	temp<-apply(apply(Tau,1,is.na), 2, sum)
    ind<-c(1:nrow(Y))[temp != 0]
  	Tau[ind,]<-rep(0,g)
    for(j in 1:length(ind))
	{
  	 Tau[ind[j], c(1:g)[logTau[ind[j],]==max(logTau[ind[j],])]]<-1
     } #j
  }
  list(Tau, logTau, reset)
} # End e-step1

