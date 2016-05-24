recluster.node.strength <-function (mat, phylo=NULL, dist="simpson", nodelab.cex=0.8, tr=100, levels=6, method="average",...){
	out<-NULL
	treesstr<-(as.phylo(hclust(recluster.dist(mat,phylo,dist),method=method)))
	tree1<-treesstr
	for (cons in 1 : levels){
		step<-0.5/(levels-1)
		p<-0.5+((cons-1)*step)
		treesstr[[cons]]<-recluster.cons(mat,phylo,dist=dist,tr=tr,method=method,p=p)$cons
	}
	btr2 <- .compressTipLabel(treesstr)
	tr2 <- recluster.check(tree1, attr(btr2, "TipLabel"))
	btr2 <- .uncompressTipLabel(btr2)
	result <- prop.clades(tr2, btr2, rooted=T)*(100/levels)
	recluster.plot(tree1,as.matrix (result), nodelab.cex = nodelab.cex,...)
	out$result<-as.matrix(result)
	out$tree<-tree1
	return(out)
}
