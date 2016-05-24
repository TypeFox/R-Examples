`chisqInd` <-
function(data,n.cat,compPval=TRUE,asMatrix=TRUE){
	comp.out<-computeContCells(data,check=FALSE,n.cat=n.cat)
	tmp<-comp.out$mat.exp
	tmp[tmp==0]<-1
	mat<-(comp.out$mat.obs-comp.out$mat.exp)^2/tmp
	out<-rowSums(mat)
	if(compPval)
		return(compChisqPval(data,out,n.cat,asMatrix=asMatrix))
	if(!asMatrix)
		return(out)
	mat<-matrix(0,nrow(data),nrow(data))
	mat[lower.tri(mat)]<-out
	if(!is.null(rownames(data)))
		colnames(mat)<-rownames(mat)<-rownames(data)
	return(mat)	
}

