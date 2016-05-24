inventorygames <-
function(n=NA,a=NA,d=NA,h=NA,m=NA,r=NA,s=NA,model=c("EOQ","EPQ")){
matriz0<-coalitions(n)[[1]]
costes<-c();costes[1]<-0
coa<-c();coa[1]<-0
for (i in 2:nrow(matriz0)){
aux<-which(matriz0[i,]!=0)
suma=0
for (j in 1:length(aux)){suma<-suma+aux[j]*10^(length(aux)-j)}
coa[i]<-suma
}
if (sum(is.na(m)!=T)==length(m)|sum(is.na(d)==T)==length(d)){ 
#caso demanda "d" desconocida
for (i in 2:nrow(matriz0)){
aux<-which(matriz0[i,]!=0)
costes[i]<-2*a*sqrt(sum(m[aux]^2))
} 
} else {
if (model=="EOQ"){
for (i in 2:nrow(matriz0)){
aux<-which(matriz0[i,]!=0)
costes[i]<-2*a*sqrt(sum(h[aux]*d[aux])/(2*a))
      } 
} 
if (model=="EPQ"){
for (i in 2:nrow(matriz0)){
aux<-which(matriz0[i,]!=0)
m<-d/sqrt(2*a*d/(h*(1-d/r))*(h+s)/s)
costes[i]<-2*a*sqrt(sum(m[aux]^2))
      } 
}

}
matriz0<-cbind(matriz0,coa,costes)
colnames(matriz0)<-c(1:n,"Coalition","Order cost")
sol<-matriz0
return(sol)}
