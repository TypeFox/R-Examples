"fun.fmkl0"<-function(k){
j<-0:k

result<-rep(NA,4)

result[1]<-0
result[2]<-(pi^2+6)/3-2
#result[2]<-integrate(function(x) (log(x)-log(1-x))^2,0,1,abs.tol=1e-100,
#stop.on.error=FALSE)$value
result[3]<-0
result[4]<-integrate(function(x) (log(x)-log(1-x))^4,0,1,abs.tol=1e-100,
stop.on.error=FALSE)$value

return(result[k])}



