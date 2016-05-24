huber.NR<-function(x,c=1.28,iter=20){

mu.k<-matrix(nrow=iter,ncol=1)
mu.k[1]<-median(x)
for(i in 1:iter){
if(mad(x) == 0)mu.k=NA
if(mad(x) != 0){
A1<-(x-mu.k[i])/mad(x)
A<-sum(sapply(A1,function(x){max(-c,min(c,x))}))
B1<-(A1>=c|A1<=-c)
B<-length(B1[B1==FALSE])
mu.k[i+1]<-mu.k[i]+((mad(x)*A)/B)
}}
mu.k
}
