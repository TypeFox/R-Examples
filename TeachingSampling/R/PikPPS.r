PikPPS<-function(n,x){
pik<- n*x/sum(x)
while((sum(pik>1))!=0){
s<-which(pik>=1)
new=(1:length(pik))[-s]
pik[s]=1
txnew<-sum(x[s])
for(k in new){
pik[k]<- (n-length(s))*x[k]/(sum(x)-txnew)
}
}
pik
}