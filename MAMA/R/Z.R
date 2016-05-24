Z<-function(lists,n.genes)
{
n.method<-length(lists)
names.method<-names(lists)
tab<-matrix(0,nrow=n.method,ncol=n.method)
colnames(tab)<-names.method
rownames(tab)<-names.method
for (i in 1:n.method)
for (j in 1:n.method)
tab[i,j]<-length(which(lists[[i]] %in% lists[[j]]))
Z<-matrix(0,nrow=n.method,ncol=n.method)
colnames(Z)<-names.method
rownames(Z)<-names.method
for (i in 1:n.method)
for (j in 1:n.method)
{Pexp<-tab[i,i]/n.genes
if (i!=j) Z[i,j]<-(tab[i,j]-tab[j,j]*Pexp)/sqrt(tab[j,j]*Pexp*(1-Pexp)) else Z[i,j]<-NA
}
return(Z)
}