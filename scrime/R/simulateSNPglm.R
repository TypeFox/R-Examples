`simulateSNPglm` <-
function(n.obs=1000,n.snp=50,list.ia=NULL,list.snp=NULL,beta0=-0.5,beta=1.5,maf=0.25,
		sample.y=TRUE,p.cutoff=0.5,err.fun=NULL,rand=NA,...){
	check.snplist<-TRUE
	if(is.null(list.snp) & is.null(list.ia)){
		list.ia<-list(c(-1,1),c(1,1,1))
		list.snp<-list(c(6,7),c(3,9,10))
		check.snplist<-FALSE
	}
	if(!is.null(list.snp) & is.null(list.ia))
		stop("If list.snp is specified, list.ia has also to be specified.")
	if(!is.list(list.ia))
		stop("list.ia must be a list.")
	n.exp<-length(unlist(list.ia))
	n.ia<-length(list.ia)
	if(is.null(list.snp)){
		list.snp<-split(1:n.exp,rep(1:n.ia,sapply(list.ia,length)))
		names(list.snp)<-names(list.ia)
		check.snplist<-FALSE
	}
	if(check.snplist){
		if(length(list.snp)!=n.ia)
			stop("list.ia and list.snp have different lengths.")
		if(any(sapply(list.ia,length)!=sapply(list.snp,length)))
			stop("Each of the vectors in list.snp must have the same length as the\n",
				"corresponding vector in list.ia.")
		if(any(!unlist(list.snp)%in%1:n.snp))
			stop("All values in list.snp must be integers between 1 and ",n.snp,".")
	} 	
	if(any(!unlist(list.ia)%in%c(-3,-2,-1,1,2,3)))
		stop("list.ia must consist of the values -3, -2, -1, 1, 2, and 3.")
	if(n.obs<10)
		stop("n.obs must be at least 10.")
	if(n.snp<max(2,n.exp))
		stop("n.snp must be at least ",max(2,n.exp),".")
	if(length(beta0)>1 | !is.numeric(beta0))
		stop("beta0 must be a numeric value.")
	if(!length(beta)%in%c(1,n.ia))
		stop("beta must either have length 1 or the same length as list.ia.")
	if(any(beta<=0))
		stop("beta must be larger than zero.")
	if(length(beta)==1)
		beta<-rep(beta,n.ia)
	if(!is.na(rand))
		set.seed(rand)
	if(!length(maf)%in%c(1,2,n.snp))
		stop("maf must be of length 1, 2 or n.snp.")
	if(any(maf<=0.01 | maf>=.5))
		stop("Any minor allele frequency must be larger than 0.01 and smaller than 0.5.")
	if(length(maf)==1)
		maf<-rep(maf,n.snp)
	if(length(maf)==2)
		maf<-runif(n.snp,maf[1],maf[2])	
	raf<-1-maf
	if(p.cutoff<=0 | p.cutoff>=1)
		stop("p.cutoff must be larger than 0 and smaller than 1.")
	mat<-matrix(0,n.obs,n.snp)
	for(i in 1:n.snp){
		tmp<-sample(0:1,2*n.obs,replace=TRUE,prob=c(raf[i],maf[i]))
		mat[,i]<-tmp[1:n.obs]+tmp[n.obs+(1:n.obs)]
	}
	mat<-mat+1
	colnames(mat)<-paste("SNP",1:n.snp,sep="")
	getIA<-function(snp,ia) 
		paste("(SNP",snp,ifelse(ia<0," != "," == "),abs(ia),")",sep="",collapse=" & ")
	vec.ia<-character(n.ia)
	for(i in 1:n.ia)
		vec.ia[i]<-getIA(list.snp[[i]],list.ia[[i]])
	mat.glm<-matrix(0,n.obs,n.ia+1)
	mat.glm[,1]<-1
	mat2<-as.data.frame(mat)
	for(i in 1:n.ia)
		mat.glm[,i+1]<-with(mat2, eval(parse(text=vec.ia[i])))
	beta<-c(beta0,beta)
	prob<-as.vector(mat.glm%*%beta)
	if(is.null(err.fun)){
		prob<-exp(prob)/(1+exp(prob))
		if(sample.y){ 
			y<-rbinom(n.obs,1,prob=prob)
			p.cutoff<-NULL
		}
		else 
			y<-ifelse(prob>p.cutoff,1,0)
		err<-NULL
		err.call<-NULL
	}
	else{
		FUN<-match.fun(err.fun)
		err<-FUN(n.obs,...)
		if(length(err)!=n.obs)
			stop("Something is wrong with err.fun.")
		y<-prob+err
		prob<-NULL
		p.cutoff<-NULL
		err.call<-getCall(match.call(expand.dots=FALSE),n.obs)
	}
	out<-list(x=mat,y=y,beta0=beta0,beta=beta[-1],ia=vec.ia,maf=maf,prob=prob,
		err=err,p.cutoff=p.cutoff,err.call=err.call)
	class(out)<-"simSNPglm"
	out			
}

