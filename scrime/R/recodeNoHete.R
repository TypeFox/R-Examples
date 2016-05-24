`recodeNoHete` <-
function(mat,geno=1:3){
	homo<-c("AA","TT","CC","GG")
	mat.homo<-matrix(0,nrow(mat),4)
	for(i in 1:4)
		mat.homo[,i]<-rowSums(mat==homo[i],na.rm=TRUE)
	rs<-rowSums(mat.homo>0)
	if(any(rs==0))
		stop("At least one of the rows of mat does neither contain homozygous\n",
			"nor heterozygous genotypes.",call.=FALSE)
	if(any(rs==1)){
		ids<-rs==1
		tmp<-mat[ids,]
		tmp[!is.na(tmp)]<-geno[1]
		mat[ids,]<-tmp
	}
	if(any(rs==2)){
		ids<-rs==2
		ref<-max.col(mat.homo[ids,,drop=FALSE])
		tmp<-mat[ids,,drop=FALSE]
		for(i in unique(ref))
			tmp[ref==i,][tmp[ref==i,]==homo[i]]<-geno[1]
		tmp[!is.na(tmp) & tmp!=1]<-geno[3]
		mat[ids,]<-tmp
	}
	mat
}

