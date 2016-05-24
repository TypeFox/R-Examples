recluster.multi <-function (tree,mat,phylo=NULL,tr=100,p=0.5,dist="simpson",method="average",boot=1000,levels=2,step=1) {
	mat<-as.matrix(mat)
	if(length(tree$tip.label)!= nrow(mat)) stop("ERROR: different site numbers between tree and matrix")

	resmulti<-matrix(nrow=nrow(mat)-1, ncol=levels)
	for(lev in 1:levels){
		resmulti[,lev]<-recluster.boot(tree,mat,phylo,tr,p,dist,method,boot,(lev-1)*step+1)
	}
	result<-resmulti
	result
}
