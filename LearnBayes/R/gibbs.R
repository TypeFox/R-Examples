gibbs=function(logpost,start,m,scale,...)
{ 
p=length(start)
vth=array(0,dim=c(m,p))
f0=logpost(start,...)
arate=array(0,dim=c(1,p))

th0=start
for (i in 1:m)
{
  for (j in 1:p)
  {
  th1=th0
  th1[j]=th0[j]+rnorm(1)*scale[j]
  f1=logpost(th1,...)
  u=runif(1)<exp(f1-f0)
  th0[j]=th1[j]*(u==1)+th0[j]*(u==0)
  f0=f1*(u==1)+f0*(u==0)
  vth[i,j]=th0[j]; 
  arate[j]=arate[j]+u
  }
}
arate=arate/m
stuff=list(par=vth,accept=arate)
return(stuff)
}


