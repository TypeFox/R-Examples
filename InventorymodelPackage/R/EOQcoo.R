EOQcoo <-
function(n=NA,a=NA,d=NA,h=NA,m=NA){
coalicion<-coalitions(n)
matriz<-as.matrix(coalicion[[1]])
matriz0<-matriz
costes<-c();costes[1]<-0
Qeoq<-matrix()
if (sum(is.na(m)!=T)==length(m)|sum(is.na(d)==T)==length(d)){ 
#caso demanda "d" desconocida
Qeoq<-EOQ(n,a,d=NA,h,m)[[1]] 
d<-Qeoq*m
}
Qeoq<-EOQ(n,a,d,h,m=NA)[[1]]
for (i in 2:nrow(matriz0)){
aux<-which(matriz0[i,]!=0)
for (j in 1:length(aux)){
matriz0[i,aux[j]]<-sqrt(2*a*d[aux[j]]^2/sum(d[aux]*h[aux]))
}
costes[i]<-a*d[aux[1]]/matriz0[i,aux[1]]+sum(h[aux]*matriz0[i,aux]/2)
}
sol<-cbind(coalicion[[2]],matriz0,costes)
rownames(sol)<-rep(" ",2^n)
print("Optimal order")
colnames(sol)<-c("Coalition",1:n,"Coalitional costs")
return(sol)}
