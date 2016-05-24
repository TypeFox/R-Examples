predictive.gpd <-
function(vector,data,threshold)
{
  data=data[data>threshold]
  linf=min(data)
  lsup=max(data)
  x=seq(linf,lsup,(lsup-linf)/70)
  n=length(x)
  int=length(vector[,1])
  res=array(0,c(n))
  for (i in 1:n)
  {for (j in 1:int)
  {res[i]=res[i]+(1/int)*dgpd(x[i],vector[j,2],threshold,vector[j,1])}}
  hist(data,freq=F,ylim=c(min(res),max(res)),breaks=seq(threshold,max(data),(max(data)-threshold)/20),main=NULL,xlab="data",ylab="density")
  lines(x,res)
  res
  
}
