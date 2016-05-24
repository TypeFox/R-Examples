`hcut` <-
function(consensus,h,group,col.text="blue",cex.text=1,...) {
dend <- consensus$dendrogram
if(class(dend)=="hclust") dend <- as.dendrogram(dend)
# cut for height
frame<-consensus$table.dend[consensus$table.dend[,4]<=h,]
cut2 <-cut(dend,h)$lower
# max clusters
max.groups<- as.numeric(length(cut2))
# cluster k
k<-group
p<-0
pclus<-rep(0,max.groups)
for(i in 1:max.groups) {
pclus[i]<-length(labels(cut2[[i]]))
}
if( k> 1) for(i in 1:(k-1))p=p+pclus[i]
#plus<-data.frame(number=pclus)
plot(cut2[[k]],...)

text(frame[,3]-p,frame[,4],labels=frame[,5],col=col.text,cex=cex.text)
return(data.frame(numbers=pclus))
}

