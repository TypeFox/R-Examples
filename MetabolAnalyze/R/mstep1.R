mstep1 <-
function(Y, Tau, Pi, mu, g)
{
  for(i in 1:g)
  {
    Pi[i]<-sum(Tau[,i])/nrow(Y)                                                                          
    mu[,i]<-apply(Y*Tau[,i],2,sum)/sum(Tau[,i])
 }
 return(list(Pi, mu))
}# End mstep1 

