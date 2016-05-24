regroup=function(data,g)
{
d=dim(data); n=d[1]; m=d[2]
N=floor(n/g)
dataG=array(0,c(N,m))
k=0
for (j in seq(1,(N-1)*g+1,g))
{
k=k+1
for (i in 0:(g-1)) 
  dataG[k,]=dataG[k,]+data[j+i,]
}
if (n>N*g)
{
for (i in (N*g+1):n)
  dataG[N,]=dataG[N,]+data[i,]
}
return(dataG)
}
