recluster.boot <- function (tree,mat,phylo=NULL,tr=100,p=0.5,dist="simpson", method="average",boot=1000,level=1) {
	mat<-as.matrix(mat)
	if(length(tree$tip.label)!= nrow(mat))
    stop("ERROR: different site numbers between tree and matrix")
	treesb<-(as.phylo(hclust(recluster.dist(mat,phylo,dist),method=method)))
	for (i in 1 : boot){
		for (testNA in 1:10000){
			xs<-mat[,sample(ncol(mat),ncol(mat)*level,replace=T)]
			if(prod(rowSums(xs))>0){
				break
			}
		}
		treesb[[i]]<-recluster.cons(xs,phylo,tr,p,dist)$cons
	}
	btr2 <- .compressTipLabel(treesb)
	tr2 <- recluster.check(tree, attr(btr2, "TipLabel"))
	btr2 <- .uncompressTipLabel(btr2)
	result <- prop.clades(tr2, btr2, rooted=T)*(100/boot)
	rows<-nrow(mat)
	matrix<-dist.nodes(tree)[(rows+1):(rows*2-1),(rows+1):(rows*2-1)]
	for (i in 2 : (rows-1)) {
		if(min(matrix[1:i-1,i])==0) {
			result[i]<- NA
		}
	}
	result<-as.matrix(result)
	result
}
