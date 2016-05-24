predictive.gumbel <-
function(vector,data)
  {
   linf=max(min(data)-1,0)
   lsup=11*max(data)/10
   x=seq(linf,lsup,(lsup-linf)/300)
   n=length(x)
   int=length(vector[,1])
   mu=vector[,1]
   sigma=vector[,2]
   res=array(0,c(n))
   for (i in 1:n)
      {for (j in 1:int)
          {
           res[i]=res[i]+(1/int)*(1/sigma[j])*exp(-(x[i]-mu[j])/sigma[j])*exp(-exp(-(x[i]-mu[j])/sigma[j]))}}
  hist(data,freq=F,ylim=c(min(res),max(res)),xlab="data",ylab="density")
  lines(x,res)
  res
  }
