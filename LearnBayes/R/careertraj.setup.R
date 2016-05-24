careertraj.setup=function(data)
{
Player=data[,1]
player.names=names(table(Player))
N=length(player.names)
m=max(table(Player))
y=array(0,c(N,m))
n=0*y
x=0*y
T=rep(0,N)
for (i in 1:N)
{
   data1=data[Player==player.names[i],]
   nk=dim(data1)
   ni=nk[1]
      for (j in 1:ni)
      {
         y[i,j]=data1[j,10]
         n[i,j]=data1[j,5]-data1[j,13] 
         x[i,j]=data1[j,3]
         T[i]=T[i]+(n[i,j]>0)
      }
}
return(list(player.names=player.names,y=y,n=n,x=x,T=T,N=N))
}