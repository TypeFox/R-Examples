Boxcox<-function(boxcox,varcod)
{
if (missing(boxcox)) 
{X2<-NULL
coup1<-NULL}
else
{
  X2<-matrix(nrow=length(varcod),ncol=length(boxcox))
  for (i in 1:length(boxcox))
  {
 	if (boxcox[i]==0)
	{
	X2[,i]<-log(varcod)}
 	else{
	X2[,i]<-(1/boxcox[i])*(varcod^{boxcox[i]}-1)}
  }
coup1<-rep(1,length(boxcox))
}
  res<-list(X2,coup1)
}