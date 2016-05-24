mct <-
function(n=NA,a=NA,av=NA,d=NA,K=NA,cooperation=c(0,1)){
if (cooperation==0){
costes<-(a+av)*d/K
sol<-costes
print("Individual cost")
}
if (cooperation==1){
matriz<-cbind(coalitions(n)[[1]],coalitions(n)[[2]],rep(0,2^n))
for (i in 2:nrow(matriz)){
aux<-which(matriz[i,1:n]==1)
matriz[i,n+2]<-max(a+av[aux])*max(d[aux]/K[aux])
}
colnames(matriz)<-c(1:n,"Coalition","Cost")
rownames(matriz)<-rep(" ",2^n)
sol<-matriz
}
return(sol)}
