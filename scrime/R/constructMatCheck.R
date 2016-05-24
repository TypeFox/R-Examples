`constructMatCheck` <-
function(n.col,vec.val,vec.equal){
	n.unequal<-sum(vec.equal==0)
	mat<-matrix(NA,2^n.unequal,n.col)
	ia.out<-ia.samp(n.unequal,2)
	geno<-0:2
	val.unequal<-vec.val[vec.equal==0]
	for(i in 1:n.unequal){
		tmp.ids<-geno[geno!=val.unequal[i]]
		ia.out[,i]<-tmp.ids[ia.out[,i]]
	}	
	mat[,vec.equal==0]<-ia.out
	ids.equal<-which(vec.equal==1)
	for(i in ids.equal)
		mat[,i]<-vec.val[i]
	mat
}

