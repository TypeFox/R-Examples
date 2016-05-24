SOC <-
function(n=NA,a=NA,d=NA,h=NA,m=NA,r=NA,s=NA,model=c("EOQ","EPQ"),cooperation=c(0,1)){
coalicion<-coalitions(n)
matriz<-as.matrix(coalicion[[1]])
matriz0<-matriz
matrizcostes<-matriz


if (sum(is.na(m)==T)!=length(m) | sum(is.na(d)!=T)==length(d)){
#Caso d conocida 
if (model=="EOQ"){Q<-EOQ(n,a,d,h,m=NA)[[1]];m=d/Q}
if (model=="EPQ"){Q<-EPQ(n,a,d,h,r,s,m=NA)[[1]];m=d/Q}
}
if (cooperation==0){coste<-2*a*m}
if (cooperation==1){
for (i in 2:nrow(matriz0)){
aux<-which(matriz0[i,]!=0)
for (j in 1:length(aux)){matrizcostes[i,aux[j]]<-2*a*m[aux[j]]^2/sqrt(sum(m[aux]^2))}
}
}

if (cooperation==0){
sol<-list(coste)
names(sol)<-c("Share the ordering costs rule (individually)")
}
if (cooperation==1){
colnames(matrizcostes)<-1:n
rownames(matrizcostes)<-rep(" ",2^n)
sol<-list(matrizcostes)

names(sol)<-c("Share the ordering costs rule (individually)")
}
return(sol)
}
