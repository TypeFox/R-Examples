covmve<-function(x){
val<-cov.mve(x)
list(center=val$center,cov=val$cov)
}
