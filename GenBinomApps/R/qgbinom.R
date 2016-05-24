qgbinom <-
function(p,size,prob,lower.tail = TRUE,log.p=FALSE){
n=sum(size)
if (log.p==TRUE){
p<-exp(p)}
if (lower.tail==FALSE){
p<-1-p}
p<-as.matrix(p)
y<-pgbinom(c(0:n),size,prob)
z<-as.matrix(apply(p,1,function(p) min(y[y>=p])))
return(apply(z,1,function(z) which(y==z,arr.ind=TRUE)[1]-1))
}
