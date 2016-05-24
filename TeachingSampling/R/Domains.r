Domains<-function(y){
y<-as.factor(y)
d<-as.double(y)
n<-length(d)
Dom<-matrix(0,n,max(d))
colnames(Dom)<-levels(y)
for(k in 1: max(d)){
Dom[,k]<-as.double(d==k)}
Dom
}