`checkMatSNP` <-
function(mat.snp,vec.val,vec.equal){
	getIdsExclude<-function(mat.snp,exclude,n.row,n.col){
		mat.check<-matrix(exclude,n.row,n.col,byrow=TRUE)
		vec.check<-rowSums(mat.snp==mat.check)
		ids<-which(vec.check==n.col)
		ids
	}
	n.row<-nrow(mat.snp)
	n.col<-length(vec.val)
	if(sum(vec.equal)==n.col)
		ids.exclude<-getIdsExclude(mat.snp,vec.val,n.row,n.col)
	else{
		ids.exclude<-NULL
		mat.exclude<-constructMatCheck(n.col,vec.val,vec.equal)
		for(i in 1:nrow(mat.exclude)){
			tmp<-getIdsExclude(mat.snp,mat.exclude[i,],n.row,n.col)
			ids.exclude<-c(ids.exclude,tmp)
		}
		ids.exclude<-unique(ids.exclude)
	}
	if(length(ids.exclude)>0)
		mat.snp<-mat.snp[-ids.exclude,,drop=FALSE]
	mat.snp
}

