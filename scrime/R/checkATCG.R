`checkATCG` <-
function(mat,first.ref=FALSE){
	homo<-c("AA","TT","CC","GG")
	hete<-c("AT","TA","AC","CA","AG","GA","TC","CT","TG","GT","CG","GC")
	if(!all(mat%in%c(homo,hete,NA)))
		stop("All values in mat must be combinations of length 2 of\n",
			"the letters A, T, C and G.")
	mat.homo<-matrix(0,nrow(mat),4)
	for(i in 1:4)
		mat.homo[,i]<-rowSums(mat==homo[i],na.rm=TRUE)
	if(any(rowSums(mat.homo>0)>2))
		stop("At least one of the SNPs shows more than 2 homozygous genotypes.",
			call.=FALSE)
	if(!first.ref){
		for(i in 2*(1:6))
			mat[mat==hete[i]]<-hete[i-1]
		hete<-hete[2*(1:6)-1]
		n.col<-6
	}
	else
		n.col<-12
	mat.hete<-matrix(0,nrow(mat),n.col)
	for(i in 1:n.col)
		mat.hete[,i]<-rowSums(mat==hete[i],na.rm=TRUE)
	if(any(rowSums(mat.hete>0)>1))
		stop("At least one of the SNPs shows more than 1 heterozygous genotype.",
			call.=FALSE)
	colnames(mat.hete)<-hete
	cs<-colSums(mat.hete)>0
	if(any(!cs))
		mat.hete<-mat.hete[,cs,drop=FALSE]
	mat.hete
}

