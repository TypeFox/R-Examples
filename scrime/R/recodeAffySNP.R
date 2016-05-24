`recodeAffySNP` <-
function(mat,refAA=FALSE,geno=1:3){
	if(!is.data.frame(mat) & !is.matrix(mat))
		stop("mat must be either a matrix or a data frame.")
	mat <- as.matrix(mat)
	if(length(geno)!=3)
		stop("geno must have length 3.")
	mat[mat=="NN"]<-NA
	if(!all(mat%in%c("AA","AB","BB",NA)))
		stop("Some of the values in mat are neither 'AA' nor 'AB' nor 'BB'.")
	mat[mat=="AB"]<-geno[2]
	if(refAA){
		mat[mat=="AA"]<-geno[1]
		mat[mat=="BB"]<-geno[3]
	}
	else{
		aa<-rowSums(mat=="AA",na.rm=TRUE)
		bb<-rowSums(mat=="BB",na.rm=TRUE)
		ids<-bb>aa
		if(sum(ids)>0){
			tmp<-mat[ids,]
			tmp[tmp=="AA"]<-geno[3]
			tmp[tmp=="BB"]<-geno[1]
			mat[ids,]<-tmp
		}
		if(sum(ids)<length(ids)){
			tmp<-mat[!ids,]
			tmp[tmp=="AA"]<-geno[1]
			tmp[tmp=="BB"]<-geno[3]
			mat[!ids,]<-tmp
		}
	}
	if(is.numeric(geno))
		mode(mat)<-"numeric"
	mat
}

