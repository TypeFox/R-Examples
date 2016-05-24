`constructMatIA` <-
function(ia.val,ia.num,ia.equal){
	n<-length(ia.val)
	mat<-matrix(NA,ia.num,n)
	geno<-0:2
	for(i in 1:n)
		mat[,i]<-if(ia.equal[i]==1) ia.val[i] 
			else sample(geno[geno!=ia.val[i]],ia.num,TRUE)
	mat
}

