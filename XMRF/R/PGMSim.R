PGMSim <-
function(n,p,alpha,Theta,maxit)
{
  X = matrix(rpois(n*p,1),n,p);
  iter = 1;
  while(iter<maxit)
    {
      par(mfrow=c(2,5))
      print(iter)
      for(j in 1:p)
        {
          X[,j] = rpois(n,exp(alpha[j] + X[,-j]%*%Theta[-j,j]))
          
        ## test and debug
        #  tt = X[,-j]%*%Theta[-j,j]
        #  sum(is.na(tt))
        #  plot(1:length(tt), tt)
        #  X[,j] = rpois(n,exp(alpha[j] + X[,-j]%*%Theta[-j,j]))
        ##tt = rpois(n,exp(alpha[j] + X[,-j]%*%Theta[-j,j]))
        }
      # sum(is.na(X))
      iter = iter + 1;
    }
  return(X)
}
