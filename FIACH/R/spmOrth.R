spmOrth<-function(a){
x<-a[,1]
for (i in 2:ncol(a)){
D = a[,i]
qr<-qr(x)
D = qr.resid(qr,D)
x<-cbind(x,D)
}
return(x)
}


