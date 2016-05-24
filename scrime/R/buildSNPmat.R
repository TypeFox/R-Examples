buildSNPmat<-function(n.obs,n.snp,maf,ias,rep=FALSE){
	n2<-2*n.obs
	raf<-1-maf
	mat<-matrix(0,n2,n.snp)
	for(i in 1:n.snp){
		tmp<-sample(0:1,2*n2,replace=TRUE,prob=c(raf[i],maf[i]))
		mat[,i]<-tmp[1:n2]+tmp[n2+(1:n2)]
	}
	mat<-mat+1
	colnames(mat)<-paste("SNP",1:n.snp,sep="")
	mat2<-as.data.frame(mat)
	mat.ias<-matrix(0,2*n.obs,length(ias))
	for(i in 1:length(ias))
		mat.ias[,i]<-with(mat2, eval(parse(text=ias[i])))
	rs<-rowSums(mat.ias)
	mat<-mat[rs<=1,]
	if(rep)
		return(mat)
	if(nrow(mat)>=n.obs){
		tmp<-sample(nrow(mat),n.obs)
		return(mat[tmp,])
	}
	k<-1
	while(nrow(mat)<n.obs){
		if(k==10){
			n.obs<-nrow(mat)
			warning("Data of only ",n.obs," observations is generated.",call.=FALSE)
			break
		}
		tmp.mat<-buildSNPmat(n.obs,n.snp,maf,ias,rep=TRUE)
		mat<-rbind(mat,tmp.mat)
		k<-k+1
	}
	tmp<-sample(nrow(mat),n.obs)
	mat[tmp,]
}


	