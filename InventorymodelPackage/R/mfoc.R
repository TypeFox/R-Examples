mfoc <-
function(n=NA,a=NA,d=NA,K=NA,cooperation=c(0,1)){
if (cooperation==0){
costs<-a*d/K
sol<-costs
print("Individual cost")
}
if (cooperation==1){
coalition<-coalitions(n)
matrix0<-as.matrix(coalition[[1]])
costs<-c();costs[1]=0
coa<-c();coa[1]=0
for (i in 2:nrow(matrix0)){
aux<-which(matrix0[i,]!=0)
sum=0
for (j in 1:length(aux)){
sum<-sum+aux[j]*10^(length(aux)-j)
}
coa[i]<-sum
costs[i]<-a*max(d[aux]/K[aux])
}
sol<-cbind(matrix0,coa,costs)
colnames(sol)<-c(1:n,"Coalition","Coalicional costs")
rownames(sol)<-rep(" ",2^n)
}
return(sol)}
