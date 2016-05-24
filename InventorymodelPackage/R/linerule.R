linerule <-
function(n=NA,a=NA,av=NA,d=NA,h=NA,m=NA){
coalicion<-STIcoo(n,a,av,d,h,m)
coalicion1<-coalicion[[1]]
coalicion2<-as.vector(coalicion1[,n+1])
matriz<-as.matrix(coalicion1[,1:n])
costes<-as.vector(coalicion1[,n+2])
matriz<-cbind(matriz,rep(0,2^n),rep(0,2^n))
for (i in 2:nrow(matriz)){
aux<-which(matriz[i,]!=0)
matriz[i,n+1]=coalicion2[i]
matriz[i,n+2]=costes[i]
}
p<-permutations(n)
aux1<-c()
for (k in 1:nrow(p)){
aux<-c()
for (i in 1:n){
for (j in 1:n){
if (i !=j){
 if(av[i]<av[j] & (which(p[k,]==i)<=which(p[k,]==j))){aux<-c(aux,1)} 
}
}
}
if (length(aux)==0){aux1<-c(aux1,k)}
}
permute<-p[aux1,]

if (is.vector(permute)==T){
phi<-c()
for (i in 1:n){
cj<-0;aux<-0
pred<-sort(permute[1:which(permute==i)-1],decreasing=F)
if (length(pred)==0) pred<-0
aux<-sort(union(pred,i),decreasing=F)
suma=0;suma1=0
for (k in 1:length(aux)){suma<-suma+aux[k]*10^(length(aux)-k)}
for (k in 1:length(pred)){suma1<-suma1+pred[k]*10^(length(pred)-k)}
cj<-cj+matriz[which(matriz[,n+1]==suma),n+2]-matriz[which(matriz[,n+1]==suma1),n+2]
phi[i]<-cj
}
}
if(is.vector(permute)==F){
phi<-c()
for (i in 1:n){
cj<-0
for (j in 1:nrow(permute)){
aux<-0
pred<-sort(permute[j,1:which(permute[j,]==i)-1],decreasing=F)
if (length(pred)==0) pred<-0
aux<-sort(union(pred,i),decreasing=F)
suma=0;suma1=0
for (k in 1:length(aux)){suma<-suma+aux[k]*10^(length(aux)-k)}
for (k in 1:length(pred)){suma1<-suma1+pred[k]*10^(length(pred)-k)}
cj<-cj+matriz[which(matriz[,n+1]==suma),n+2]-matriz[which(matriz[,n+1]==suma1),n+2]
}
phi[i]<-cj/(nrow(permute))
}
}
return(phi)}
