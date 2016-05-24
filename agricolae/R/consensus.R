`consensus` <-
function(data,distance=c("binary","euclidean","maximum","manhattan","canberra",
"minkowski","gower"),method=c("complete","ward","single","average","mcquitty","median",
"centroid"),nboot=500,duplicate=TRUE,cex.text=1,col.text="red", ...)
{
t0<-Sys.time()
distance <- match.arg(distance)
method <- match.arg(method)
if(distance=="gower") {
if (requireNamespace("cluster", quietly = TRUE)) {
distancia<-cluster::daisy(data,metric=distance)
}
else{
return("Please install cluster package for calculate gower distance")
}
}
else distancia<-dist(data,method=distance)
nc<-ncol(data)
nr<-nrow(data)
dend<-hclust(distancia,method=method)
h1<-cutree(dend,h=0)
h2<-data.frame(h1)
h3<-unique(h2)
dup<-duplicate
duplicate<-NULL
# To study duplicate
if(dup){
if(nrow(h3) < length(h1)){
nr<-nrow(h3)
data<-data.frame(d=rownames(data),data)
h3<-data.frame(d=rownames(h3),h3)
duplicate<- merge(h3,data,by="d",all=TRUE)
duplicate<-duplicate[is.na(duplicate[,2]),]
dup0<-duplicate[,-2]
duplicate<-as.character(duplicate$d)
data<-merge(h3,data,by="d")
dup1 <-data[,-2]
dup0<-cbind(dup0,unique="")
dup0[,1]<-as.character(dup0[,1])
dup1[,1]<-as.character(dup1[,1])
ncdup<-ncol(dup1)
dup0[,ncdup+1]<-as.character(dup0[,ncdup+1])
ndup0<-nrow(dup0)
ndup1<-nrow(dup1)
ncdup<-ncol(dup1)
for ( i in 1:ndup0) {
for ( j in 1:ndup1) {
if(sum(dup0[i,2:ncdup]==dup1[j,-1],na.rm=TRUE) == ncdup-1){
dup0[i,ncdup+1]<-dup1[j,1]
break
}
}
}
if (sum(dup0[,ncdup+1]=="")>0) {
add1<-dup0[dup0[,ncdup+1]=="",]
add1<-data.frame(d=add1[,1],h1=0,add1[,2:ncdup])
data<-rbind(data,add1)
}
rownames(data)<-data[,1]
data<-data[,c(-1,-2)]
nc<-ncol(data)
if(distance=="gower") distancia<-daisy(data,metric=distance)
else distancia<-dist(data,method=distance)
dend<-hclust(distancia,method=method)
duplicate<-dup0[dup0[,ncdup+1]!="",][,c(1,ncdup+1)]
names(duplicate)[1]<-"duplicate"
}}
nr<-nrow(data)
if(!is.null(duplicate)) {
cat("\nDuplicates:",nrow(duplicate))
cat("\nNew data  :", nr,"Records\n")
}

dinicio<-dend$merge
d0<-hgroups(dinicio)
clases<- data.frame(d=d0,height=dend$height,sq=1:length(d0))
#rownames(clases)<-d0

b<- nboot
d<-NULL
#########################
for ( k in 1: b) {
muestra<-sample(1:nc,replace=TRUE)
boot1<-data[,muestra]
if(distance=="gower") distancia<-daisy(boot1,metric=distance)
else distancia<-dist(boot1,method=distance)
d1<-hclust(distancia,method=method)$merge
d1<-hgroups(d1)
d<-rbind(d,d1)
}
td<-table(d)
td<-data.frame(td)
junto<-merge(clases,td,by="d")
junto[,4]<-junto[,4]*100/b
junto<-merge(clases,junto,by="d",all=TRUE)
junto[is.na(junto)]<-0
junto<-junto[order(junto[,3]),]
############
tiempo<- Sys.time()-t0
unidad <- attributes(tiempo)$units
cat("\nConsensus hclust\n" )
cat("\nMethod distance:",distance)
cat("\nMethod cluster :",method)
cat("\nrows and cols  :",nr,nc)
cat("\nn-bootstrap    :",nboot)
cat("\nRun time       :",tiempo,unidad,"\n\n")
cc<-cbind(dend$merge,height=dend$height,percentage=round(junto[,6],1))
co<-dend$order
n1<-nrow(cc)
n2<-length(co)
p<-rep(0,n1)
for(i in 1:n1) {
if ((cc[i,1] < 0) & (cc[i,2] < 0) ) {
for(k in 1:n2) {
if(co[k]==-cc[i,1]) k1=k
if(co[k]==-cc[i,2]) k2=k
}
p[i]<- (k1+k2)/2
}
if ((cc[i,1]) < 0 & (cc[i,2] > 0)) {
for(k in 1:n2) {
if(co[k]==-cc[i,1]) k1=k
}
p[i]<-(k1+p[cc[i,2]])/2
}
if ((cc[i,1] > 0) & (cc[i,2] < 0)) {
for(k in 1:n2) {
if(co[k]==-cc[i,2]) k1=k
}
p[i]<-(k1+p[cc[i,1]])/2
}
if ((cc[i,1] > 0) & (cc[i,2] > 0)) {
p[i]<-(p[cc[i,1]]+p[cc[i,2]])/2
}
}
table.dend <- data.frame(dend$merge,xaxis=p,cc[,3:4],groups=d0)
plot(dend,...)
#delta<-0.05*(max(cc[,4]))
text(p,cc[,3],ceiling(cc[,4]),cex=cex.text,col=col.text)
return(list(table.dend=table.dend,dendrogram=dend,duplicates=duplicate) )
}
