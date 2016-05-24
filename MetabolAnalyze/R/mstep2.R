mstep2 <-
function(Y, Tau, Pi, mu, W, Sig, g, p, q)
{
  N<-nrow(Y)
  Vp<-0.1
  C2p<-0.1
    
  sig<-rep(NA,g)
  for(i in 1:g)
  {
     yc<-sweep(Y,2,mu[,i],"-")
     M_1<-solve((t(W[,,i])%*%W[,,i]) + (Sig)*diag(q))
     SW<-(t(yc*Tau[,i])%*%(yc%*%W[,,i]))/sum(Tau[,i])

     W[,,i]<-SW%*%solve((Sig)*diag(q) + ((M_1%*%t(W[,,i]))%*%SW))                                            
     sig[i]<-Pi[i]*(sum(diag(((t(yc*Tau[,i])%*%yc)/sum(Tau[,i])) - (SW%*%(M_1%*%t(W[,,i]))))))               
  } 
  MLESig<-(1/p)*sum(sig)
  
  Sig<- ((MLESig*p) + C2p)*(N/(Vp + 2 + (N*p)))
  return(list(W, Sig))
}# End mstep2

