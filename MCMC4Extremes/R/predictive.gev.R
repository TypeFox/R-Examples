predictive.gev <-
function(vector,data)
{
  linf=max(min(data)-1,0)
  lsup=11*max(data)/10
  x=seq(linf,lsup,(lsup-linf)/70)
  n=length(x)
  int=length(vector[,1])
  res=array(0,c(n))
  for (i in 1:n)
  {for (j in 1:int)
  {
    
    if ((vector[j,3]>0) && (x[i]>(vector[j,1]-vector[j,2]/vector[j,3])))      
      res[i]=res[i]+(1/int)*dgev(x[i],vector[j,3],vector[j,1],vector[j,2])
    if ((vector[j,3]<0) && (x[i]<(vector[j,1]-vector[j,2]/vector[j,3]))) 
      res[i]=res[i]+(1/int)*dgev(x[i],vector[j,3],vector[j,1],vector[j,2])   
  }}
  hist(data,freq=F,ylim=c(min(res),max(res)),main=NULL,xlab="data",ylab="density")
  lines(x,res)
  res
  
}
