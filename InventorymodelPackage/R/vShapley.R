vShapley <-
function(n,game){
costes<-characteristicfunction(n,game)
#library(e1071) #requiere paquete e1071
p<-permutations(n)
phi<-c()
for (i in 1:n){
cj<-0
for (j in 1:nrow(p)){
aux<-0
pred<-sort(p[j,1:which(p[j,]==i)-1],decreasing=F)
if (length(pred)==0) pred<-0
aux<-sort(union(pred,i),decreasing=F)
suma=0;suma1=0
for (k in 1:length(aux)){suma<-suma+aux[k]*10^(length(aux)-k)}
for (k in 1:length(pred)){suma1<-suma1+pred[k]*10^(length(pred)-k)}
cj<-cj+costes[which(costes[,n+1]==suma),n+2]-costes[which(costes[,n+1]==suma1),n+2]
}
phi[i]<-cj/nrow(p)
}
return(phi)}
