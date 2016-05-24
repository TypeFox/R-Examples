probs.to.matrix<-function(probs, genenames)
{
xx<-matrix(FALSE, nrow=length(genenames), ncol=length(names(probs)))
for (i in 1:length(probs))
xx[,i]<-genenames %in% probs[[i]]
colnames(xx)<-names(probs)
rownames(xx)<-genenames
return(xx)
}
