recluster.cons <- function(mat,phylo=NULL,tr=100,p=0.5,dist="simpson", method="average", blenghts=TRUE, select=FALSE) {
      if(data.class(mat)=="dist")
      {distance<-mat}
      else
      {distance<-recluster.dist(mat,phylo,dist)}
      sampl<-cbind((1:nrow(as.matrix(distance)))+10,rownames(as.matrix(distance)))
      res<-NULL	
      trees<-NULL
      RSS<-NULL
      for (i in 1 : tr){
		dist1<-as.matrix(distance)		
		sampl[,1]<-sample(1:nrow(sampl))+10
		rownames(dist1)<-sampl[,1]
		colnames(dist1)<-sampl[,1]
		dist1<-dist1[order(rownames(dist1)), order(colnames(dist1))]
		nam<-sampl[order(sampl[,1]),]
		rownames(dist1)<-nam[,2]
		tree<-as.phylo(hclust(as.dist(dist1),method=method))
		if (select){RSS[[i]]<-attr(nnls.tree(as.dist(dist1), tree, rooted=T),"RSS")}
		trees[[i]]<-tree				
	}
	if (select){trees[which(RSS<median(RSS))]}
	cons<-compute.brlen(consensus(trees[1:i],p=p, check.labels=T), method="Grafen")
	if(blenghts){cons<-nnls.tree(distance,cons, rooted=T,trace=F)}
	res$cons<-multi2di(cons, random=T)
	res$trees<-trees
	res$RSS<-RSS
	return(res)
}
