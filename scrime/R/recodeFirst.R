`recodeFirst` <-
function(mat,bases,geno=1:3){
	oldgeno<-paste(rep(bases,e=2),rep(bases,2),sep="")[-3]
	for(i in 1:3)
		mat[mat==oldgeno[i]]<-geno[i]
	mat
}

