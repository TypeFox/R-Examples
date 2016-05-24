`recodeCount` <-
function(mat,bases,geno=1:3){
	oldgeno<-paste(rep(bases,e=2),rep(bases,2),sep="")
	aa<-rowSums(mat==oldgeno[1],na.rm=TRUE)
	bb<-rowSums(mat==oldgeno[4],na.rm=TRUE)
	ids<-bb>aa
	if(sum(ids)>0){
		tmp<-mat[ids,]
		tmp[tmp==oldgeno[4]]<-geno[1]
		tmp[tmp==oldgeno[1]]<-geno[3]
		mat[ids,]<-tmp
	}
	if(sum(ids)<length(ids)){
		tmp<-mat[!ids,]
		tmp[tmp==oldgeno[4]]<-geno[3]
		tmp[tmp==oldgeno[1]]<-geno[1]
		mat[!ids,]<-tmp
	}
	mat[mat%in%oldgeno[2:3]]<-geno[2]
	mat
}

