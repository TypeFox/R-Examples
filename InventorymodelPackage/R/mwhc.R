mwhc <-
function(n=NA,a=NA,b=NA,d=NA,K=NA,cooperation=c(0,1),allocation=c(0,1)){
if (sum(sort(d/K)==(d/K))!=n) {
print("Warning: agents, d and K are not in the order indicated by the ratios d/K.")
d0K0<-d/K
d<-d[order(d0K0)]
K<-K[order(d0K0)]
}
if (cooperation==0){
pedido<-sqrt(b*d/(2*a+b*K^2/d))
costes<-sqrt(b*d*(2*a+b*K^2/d))-b*K
faltantes<-d/pedido-K
sol<-list(pedido,costes,faltantes)
names(sol)<-c("Number of orders per time unit","Costs","Optimal shortages")
}
if (cooperation==1){
matriz0<-coalitions(n)
matriz<-cbind(as.vector(matriz0[[2]]),rep(0,2^n),rep(0,2^n))
matriz0<-matriz0[[1]]
matriz1<-matrix(0,ncol=n+1,nrow=2^n)
for (i in 2:nrow(matriz)){
matriz0[i,]->coa
aux<-which(coa==1)
coaux<-which(coa==1)
s=length(coaux)
k=s+1
T<-rep(0,length(coa))
xT=0
SxT=coa

while(sum(SxT==T)!=length(coa)){
k=k-1;k
T[coaux[k]]=1;T
aux1<-which(T==1)
xT<-sqrt(sum(b[aux1]*d[aux1])/(2*a+sum(b[aux1]*K[aux1]^2/d[aux1])))

aux2<-which(xT<d/K)
SxT<-rep(0,length(coa));SxT[intersect(which(coa==1),aux2)]=1;SxT
}
matriz[i,2]<-xT
ind<-which(SxT==1)
matriz[i,3]<-sum(b[ind]*d[ind])/xT-sum(b[ind]*K[ind])
suma1<-0
c0<-sqrt((2*a+sum(b[ind]*K[ind]^2/d[ind]))/sum(b[ind]*d[ind]))
if (allocation==1){
for (k in 1:length(ind)){
suma1<-suma1+ind[k]*10^(length(ind)-k)
matriz1[i,1+ind[k]]<-c0*b[ind[k]]*d[ind[k]]-b[ind[k]]*K[ind[k]]
}
matriz1[i,1]<-suma1
}
}
colnames(matriz)<-c("Coalitions","Optimal orders","Costs")
sol<-matriz
rownames(matriz)<-rep(" ",2^n)
colnames(matriz1)<-c("Coalition_SxT",1:n)
rownames(matriz1)<-rep(" ",2^n)
if (allocation==1) {
sol<-list(matriz,matriz1);names(sol)<-c("Optimal policies","R-rule")
}
}
return(sol)}
