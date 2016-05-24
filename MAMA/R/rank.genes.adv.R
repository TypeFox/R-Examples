rank.genes.adv<-function(testp)
{
n<-ncol(testp$test)
for (i in 1:n)
{
 if (i==1) {
 r<-rank.genes(testp$test[,i],testp$p[,i])
 } else {
 r<-cbind(r,rank.genes(testp$test[,i],testp$p[,i])) 
 }
}
rownames(r)<-rownames(testp$test)
colnames(r)<-colnames(testp$p)
return(r)
}