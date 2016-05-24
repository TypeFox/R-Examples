`recodeSNPs` <-
function(mat,first.ref=FALSE,geno=1:3,snp.in.col=FALSE){
	if(!is.data.frame(mat) & !is.matrix(mat))
		stop("mat must be either a matrix or a data frame.")
	mat<-as.matrix(mat)
	if(length(geno)!=3)
		stop("geno must have length 3.")
	if(snp.in.col)
		mat<-t(mat)
	mat[mat=="NN"]<-NA
	mat.hete<-checkATCG(mat,first.ref=first.ref)
	ids.nohete<-rowSums(mat.hete)==0
	if(any(!ids.nohete)){
		cn<-strsplit(colnames(mat.hete),"")
		FUN<-if(first.ref) recodeFirst else recodeCount
		for(i in 1:length(cn)){
			ids<-mat.hete[,i]>0
			mat[ids,]<-FUN(mat[ids,,drop=FALSE],cn[[i]],geno=geno) 
		}
	}
	if(any(ids.nohete))
		mat[ids.nohete,]<-recodeNoHete(mat[ids.nohete,, drop=FALSE],geno=geno)
	if(is.numeric(geno))
		mode(mat)<-"numeric"
	if(snp.in.col)
		mat<-t(mat)
	mat
}

