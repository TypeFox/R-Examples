pituus<-function(x){
#laskee euklid pituuden nelion matriisien x riveille
#
d<-length(x[1,])
lkm<-length(x[,1])
vast<-matrix(0,lkm,1)
i<-1
while (i<=lkm){
  j<-1
  while (j<=d){
    vast[i]<-vast[i]+(x[i,j])^2
    j<-j+1
  }
  i<-i+1
}
return(t(vast))
}
