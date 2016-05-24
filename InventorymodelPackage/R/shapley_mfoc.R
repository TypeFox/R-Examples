shapley_mfoc <-
function(n=NA,a=NA,d=NA,K=NA){
dk<-order(d/K)
d<-d[dk];K<-K[dk]
cind<-as.vector(mfoc(n,a,d,K,cooperation=0))
shapley<-c();shapley[1]<-cind[1]/n
for (i in 2:n){
aux<-0
for (j in 2:i){aux<-aux+(cind[j]-cind[j-1])/(n-j+1)}
shapley[i]<-shapley[1]+aux
}
print("Shapley value")
return(shapley)}
