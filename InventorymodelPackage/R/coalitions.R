coalitions <-
function(n){
coa<-matrix(0,nrow=1,ncol=n)
while(nrow(coa)<2^n){
aux<-sample(c(0,1),n,replace=T)
dis=0
for (i in 1:nrow(coa)){
   if(sum(aux==coa[i,])==n){break} else {dis=dis+1}
   if (dis==nrow(coa)){coa<-rbind(coa,aux)}
   }
}
aux1<-c();aux1[1]<-0
for (i in 2:nrow(coa)){
aux<-which(coa[i,]!=0)
sum=0
for (j in 1:length(aux)){
sum<-sum+aux[j]*10^(length(aux)-j)
}
aux1<-c(aux1,sum)
}
coa<-coa[order(aux1),]

coalitions<-sort(aux1)
rownames(coa)<-c()
colnames(coa)<-c()
sol<-list(coa,coalitions)
names(sol)<-c("Binary","Classic")
return(sol)}
