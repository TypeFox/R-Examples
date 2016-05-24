EPQcoo <-
function(n=NA,a=NA,d=NA,h=NA,m=NA,r=NA,s=NA){
coalicion<-coalitions(n)
matriz<-as.matrix(coalicion[[1]])
matriz0<-matriz
matrizf<-matriz
costes<-c();costes[1]<-0
if (sum(is.na(m)==T)!=length(m)|sum(is.na(d)==T)==length(d)){ 
#caso demanda "d" desconocida
Qepq<-EPQ(n,a,d=NA,h,m,r,s)[[1]]
d<-Qepq*m
} #else {
Qepq<-EPQ(n,a,d,h,m=NA,r,s)[[1]]
for (i in 2:nrow(matriz0)){
aux<-which(matriz0[i,]!=0)
for (j in 1:length(aux)){
matriz0[i,aux[j]]<-sqrt(2*a*d[aux[j]]^2/sum(d[aux]*h[aux]*s[aux]*(1-d[aux]/r[aux])/(h[aux]+s[aux])))
matrizf[i,aux[j]]<-matriz0[i,aux[j]]*h[aux[j]]*(1-d[aux[j]]/r[aux[j]])/(h[aux[j]]+s[aux[j]])
}
costes[i]<-2*a*sqrt(sum((d[aux]/Qepq[aux])^2))
}
#}
matriz0<-cbind(matriz0,costes)
rownames(matriz0)<-rep("",2^n)
rownames(matrizf)<-rep("",2^n)
colnames(matriz0)<-c(1:n,"Costs")
colnames(matrizf)<-1:n
sol<-list(matriz0,matrizf)
names(sol)<-c("Optimal order","Optimal shortages")

return(sol)}
