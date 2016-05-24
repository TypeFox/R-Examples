ositalog<-function(n,osimat,sara){
#
rownum<-dim(osimat)[1]
#rividim<-dim(osimat)[2]
#
#Lasketaan otoskoko
#las<-0
#j<-1
#while (j<=rividim)
#  if (osimat[saradim,j]!=NaN) las<-las+1
#  endif
#  j<-j+1
#endo
#n<-(saradim-1)*rividim+las
#Ryhdytaan paatehtavaan
#
tulos<-matrix(0,n,1)
fal<-0
tru<-1
#
#i=1
#while (i<=n)
#   tulos(i,1)<-false
#endo
#
i<-1
while (i<=rownum){
  if (!is.na(osimat[i,sara])) tulos[osimat[i,sara]]<-tru 
  i<-i+1
}
return(tulos)
}
