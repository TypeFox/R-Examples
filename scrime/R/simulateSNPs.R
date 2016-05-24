`simulateSNPs` <-
function(n.obs,n.snp,vec.ia,prop.explain=1,list.ia.val=NULL,vec.ia.num=NULL,vec.cat=NULL,
		maf=c(0.1,0.4),prob.val=rep(1/3,3),list.equal=NULL,prob.equal=0.8,
		rm.redundancy=TRUE,shuffle=FALSE,shuffle.obs=FALSE,rand=NA){
	raf<-NULL
	if(!length(n.obs)%in%1:2)
		stop("n.obs must be either of length 1 or 2.")
	if(!is.null(vec.cat)){
		if(length(vec.cat) != length(vec.ia))
			stop("vec.cat must have the same length as vec.ia.")
		n.cat <- length(unique(vec.cat)) + 1
	}
	else
		n.cat <- 2
	if(length(n.obs) == 1){
		n.control <- ceiling(n.obs / n.cat)
		n.case <- (n.cat - 1) * n.control
	}
	else{
		n.case<-n.obs[1]
		n.control<-n.obs[2]
	}
	n.obs<-n.case+n.control
	if(n.case<10 | n.control<10)
		stop("Both the number of cases and the number of controls should be at least 10.")
	if(!is.na(rand))
		set.seed(rand)
	if(any(vec.ia<1))
		stop("Each entry of vec.ia must be at least 1.")
	n.ia<-length(vec.ia)
	if(n.snp<sum(vec.ia))
		stop("sum(vec.ia) must be <= n.snp.")
	if(!length(prop.explain)%in%c(1,n.ia))
		stop("prop.explain must be either of length 1 or length(vec.ia).")
	if(any(prop.explain>1 | prop.explain<=.5))
		stop("prop.explain must be >0.5 and <=1.")
	if(length(prop.explain)==1)
		prop.explain<-rep(prop.explain,n.ia)
	if(is.null(list.ia.val)){
		list.ia.val<-vector("list",n.ia)
		for(i in 1:n.ia)
			list.ia.val[[i]]<-sample(0:2,vec.ia[i],replace=TRUE,prob=prob.val)
	}
	if(length(list.ia.val)!=n.ia)
		stop("list.ia.val must have the same length as vec.ia.")
	tmp<-unlist(lapply(list.ia.val,length))
	if(any(tmp!=vec.ia))
		stop("The length of the vectors in list.ia.val must correspond to the orders of",
			" the interactions in vec.ia.")
	if(any(!unlist(list.ia.val)%in%0:2))
		stop("All entries in list.ia.val must be either 0, 1 or 2.")
	if(prob.equal<0 | prob.equal>1)
		stop("prob.equal must be between 0 and 1.")
	if(is.null(list.equal)){
		list.equal<-vector("list",n.ia)
		for(i in 1:n.ia)
			list.equal[[i]]<-sample(0:1,vec.ia[i],replace=TRUE,prob=c(1-prob.equal,
				prob.equal))
	}
	if(length(list.equal)!=n.ia)
		stop("list.equal must have the same length as vec.ia.")
	tmp<-unlist(lapply(list.equal,length))
	if(any(tmp!=vec.ia))
		stop("The length of the vectors in list.equal must correspond to the orders of",
			" the interactions specified in vec.ia.")
	if(any(!unlist(list.equal)%in%0:1))
		stop("All entries in list.equal must be either 0 or 1.")
	if(is.null(vec.ia.num)){
		vec.ia.num<-rep(n.case%/%n.ia,n.ia)
		tmp<-sample(n.ia,n.case%%n.ia)
		vec.ia.num[tmp]<-vec.ia.num[tmp]+1
	}
	if(length(vec.ia.num)!=n.ia)
		stop("vec.ia.num must have the same length as vec.ia.")
	if(sum(vec.ia.num)>n.case)
		stop("sum(vec.ia.num) is larger than the number of cases.")
	if(any(vec.ia.num<10))
		stop("All the entries in vec.ia.num must be >= 10.")
	if(!length(maf)%in%c(1,2,n.snp))
		stop("maf must be of length 1, 2 or n.snp.")
	if(any(maf<=0.01 | maf>=.5))
		stop("Any minor allele frequency must be larger than 0.01 and smaller than 0.5.")
	if(length(maf)==1)
		maf<-rep(maf,n.snp)
	if(length(maf)==2)
		maf<-runif(n.snp,maf[1],maf[2])	
	if(is.null(raf))
		raf<-1-maf
	else{
		if(!length(raf)%in%c(1,n.snp))
			stop("raf must be either of length 1 or n.snp.")
		if(any(raf<=0.5 | raf>=0.99))
			stop("Any entry of raf must be larger than 0.5 and smaller than 0.99.")
		if(length(raf)==1)
			raf<-rep(raf,n.snp)
	}
	vec.ia.control<-vec.ia.num/prop.explain*(1-prop.explain)
	vec.ia.control<-round(vec.ia.control)
	if(sum(vec.ia.control)>n.control)
		stop("The number of controls explained by the interactions is larger than the",
			" total number of controls.")
	mat<-matrix(NA,n.obs,n.snp)
	vec.ia.cum<-c(0,cumsum(vec.ia))
	vec.case.cum<-c(0,cumsum(vec.ia.num))
	vec.control.cum<-c(0,cumsum(vec.ia.control))+n.case
	for(i in 1:n.ia){
		tmp.obs<-tmp.case<-(vec.case.cum[i]+1):vec.case.cum[i+1]
		tmp.ia<-(vec.ia.cum[i]+1):vec.ia.cum[i+1]
		mat[tmp.case,tmp.ia]<-constructMatIA(list.ia.val[[i]],vec.ia.num[i],list.equal[[i]])
		if(vec.ia.control[i]!=0){
			tmp.control<-(vec.control.cum[i]+1):vec.control.cum[i+1]
			mat[tmp.control,tmp.ia]<-constructMatIA(list.ia.val[[i]],vec.ia.control[i],
				list.equal[[i]])
			tmp.obs<-c(tmp.obs,tmp.control)
		}
		n.not.exp<-n.obs-length(tmp.obs)
		mat[-tmp.obs,tmp.ia]<-sample.snpmat(n.not.exp,vec.ia[i],maf[tmp.ia],raf[tmp.ia],
			list.ia.val[[i]],list.equal[[i]])
	}
	n.iasnp<-sum(vec.ia)
	if(n.snp>n.iasnp){
		for(i in (n.iasnp+1):n.snp){
			tmp<-sample(0:1,2*n.obs,replace=TRUE,prob=c(raf[i],maf[i]))
			mat[,i]<-tmp[1:n.obs]+tmp[(n.obs+1):(2*n.obs)]
		}
	}
	y<-rep(c(1,0),c(n.case,n.control))
	if(shuffle.obs){
		ids.obs<-sample(n.obs)
		y<-y[ids.obs]
		mat<-mat[ids.obs,]
	}
	if(shuffle){
		ids.snp<-sample(n.snp)
		mat<-mat[,ids.snp]
	}
	else
		ids.snp<-1:n.snp
	colnames(mat)<-paste("SNP",1:n.snp,sep="")
	rownames(mat)<-1:n.obs
	ias<-paste("SNP",ids.snp[1:n.iasnp],sep="")
	tmp.sign<-ifelse(unlist(list.equal)==0,"!=","==")
	ias<-paste(ias,tmp.sign,unlist(list.ia.val))
	andor<-rep(rep(c("  &  ","  OR  "),n.ia),as.numeric(rbind(vec.ia-1,1)))
	ias<-paste(ias,andor,sep="",collapse="")
	ias<-unlist(strsplit(ias,"  OR  "))
	ias<-removeRedundancy(mat,ias)
	mat2<-as.data.frame(mat)
	vec.cases<-vec.controls<-numeric(n.ia)
	mat.eval <- matrix(0, n.obs, n.ia)
	for(i in 1:n.ia)
		mat.eval[,i] <- with(mat2, eval(parse(text=ias[i])))
	vec.cases <- colSums(mat.eval[y==1,])
	vec.controls <- colSums(mat.eval[y==0,])
	if(any(vec.cases!=vec.ia.num))
		stop("Something went wrong when removing redundancy SNPs.")
	if(any(vec.controls!=vec.ia.control))
		stop("Something went wrong when removing redundancy SNPs.")
	if(is.null(vec.cat))
		tab.explained <- data.frame(Interaction=ias, Cases=vec.cases,
			Controls=vec.controls, stringsAsFactors=FALSE)
	else{
		mat.eval <- mat.eval[y==1,]
		idsSeveral <- rowSums(mat.eval) > 1
		if(any(idsSeveral))
			stop("This should not happen. Please inform the author of this function.")
		y1 <- max.col(mat.eval)
		y[y==1] <- vec.cat[y1]		 
		tab.explained <- data.frame(Interaction=ias, Level=vec.cat, Cases=vec.cases,
			Controls=vec.controls, stringsAsFactors=FALSE)
	}
	out<-list(data=mat,cl=y,tab.explain=tab.explained,ia=ias,maf=maf)
	class(out)<-"simulatedSNPs"
	out	
}

