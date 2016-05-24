taulc<-function(x,mu.too=F){
#
val<-tauvar(x)
if(mu.too){
val[2]<-val
val[1]<-tauloc(x)
}
val
}

