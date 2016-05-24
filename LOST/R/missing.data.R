missing.data <-
function(x,remperc) {
x<-as.matrix(x)
totaldata<-nrow(x)*ncol(x)
n<-round(totaldata*remperc)
ndat<-1:totaldata
remove<-sample(ndat,n,replace=FALSE)
for (k in 1:n) {
  i<-remove[k]
	x[i]<-NA
}
return(x)
}
