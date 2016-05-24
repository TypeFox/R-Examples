`Get1G` <-
function(phi,n){
 p<-length(phi)
 x0<-sum(c(1,-phi))^2
 x<-x0-rowSums(GetB(phi))-GetKappa(phi)
 c(x,rep(x0,n-2*p),rev(x))
}

