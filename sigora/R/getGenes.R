getGenes<-function(yy,i){
idmap<-get(data(idmap,envir=as.environment(parent.frame())))
yy1<-yy$detailes_results[which(yy$detailes_results[,'pathway']==yy$summary_results[i,1]),,drop=FALSE]
## pgs<-unique(as.character(as.vector(yy1[yy1[,4]==1,1:2])))
## print(pgs)
r1<-aggregate(as.numeric(yy1[,4]),by=list(yy1[,1]),FUN=sum)
r2<-aggregate(as.numeric(yy1[,4]),by=list(yy1[,2]),FUN=sum)
r3<-aggregate(rbind(r1,r2)[,2],by=list(rbind(r1,r2)[,1]),FUN=sum)
r4<-as.data.frame(r3[order(-r3[,2]),])
r4[,2]<-0.5*r4[,2]
colnames(r4)<-c("gene","contribution")
v1<-(sapply(idmap,function(x)length(intersect(x,r4[,1]))))
if(max(v1)>0){
t1<-which.max(v1)
r4<-cbind(r4,idmap[match(r4[,1],idmap[,t1]),-t1])
}        
r4
}
