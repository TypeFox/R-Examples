make.matrix<-function(lists)
{
n.method<-length(lists)
names.method<-names(lists)
temp1<-lists[[1]]
for (i in 2:n.method) {temp1<-c(temp1,lists[[i]])}
genes<-unique(temp1)

HeatMatr<-c(rep(0,length(genes)))
HeatMatr[genes %in% lists[[1]]]<-1
for (i in 2:n.method)
{
temp3<-c(rep(0,length(genes)))
temp3[genes %in% lists[[i]]]<-1
HeatMatr<-cbind(HeatMatr,temp3)
}
rownames(HeatMatr)<-genes
colnames(HeatMatr)<-names.method
return(HeatMatr)
}
