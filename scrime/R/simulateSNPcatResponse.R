simulateSNPcatResponse<-function(n.obs=1000,n.snp=50,list.ia=NULL,list.snp=NULL,withRef=FALSE,
		beta0=-0.5,beta=1.5,maf=0.25,sample.y=TRUE,rand=NA){
	check.snplist<-TRUE
	if(is.null(list.snp) & is.null(list.ia)){
		list.ia<-list(c(-1,1),c(1,1,1),list(c(-1,1),c(1,1,1)))
		list.snp<-list(c(6,7),c(3,9,10),list(c(2,5),c(1,4,8)))
		check.snplist<-FALSE
	}
	if(!is.null(list.snp) & is.null(list.ia))
		stop("If list.snp is specified, list.ia has also to be specified.")
	if(!is.list(list.ia))
		stop("list.ia must be a list.")
	n.exp<-length(unlist(list.ia))
	n.cat<-length(list.ia)
	if(n.cat<2)
		stop("list.ia must consist of at least two objects.")
	n.ias<-sapply(list.ia,function(x) if(is.list(x)) length(x) else 1)
	n.vars<-sapply(list.ia,function(x) length(unlist(x)))
	if(is.null(list.snp)){
		list.snp<-split(1:n.exp,rep(1:n.cat,n.vars))
		names(list.snp)<-names(list.ia)
		if(any(n.ias>1)){
			tmp.ids<-which(n.ias>1)
			for(i in tmp.ids)
				list.snp[[i]]<-split(list.snp[[i]],rep(1:n.ias[i],
					sapply(list.ia[[i]],length)))
		}
		check.snplist<-FALSE
	}
	if(check.snplist){
		if(length(list.snp)!=n.cat)
			stop("list.ia and list.snp have different lengths.")
		if(any(!unlist(list.snp)%in%1:n.snp))
			stop("All values in list.snp must be integers between 1 and ",n.snp,".")
		tmp<-sapply(list.snp,function(x) length(unlist(x)))
		if(any(tmp!=n.vars))
			stop("The lengths of the objects in list.ia and list.snp differ.")
		if(any(n.ias>1)){
			tmp.ids<-which(n.ias>1)
			for(i in tmp.ids){
				if(any(sapply(list.ia[[i]],length)!=sapply(list.snp[[i]],length)))
					stop("list.ia and list.snp must have the same structure.")
			}
		}
	}
	if(any(!unlist(list.ia)%in%c(-3,-2,-1,1,2,3)))
		stop("list.ia must consist of the values -3, -2, -1, 1, 2, and 3.")
	if(n.obs<20)
		stop("n.obs must be at least 20.")
	if(n.snp<max(2,n.exp))
		stop("n.snp must be at least ",max(2,n.exp),".")
	if(length(beta0)==1)
		beta0<-rep(beta0,n.cat)
	if(length(beta0)!=n.cat)
		stop("beta0 must be either a numeric value or a vector of the same length as list.ia.")
	if(any(unlist(beta)<=0))
		stop("All values in beta must be larger than 0.")
	if(length(beta)==1)
		beta<-split(rep(beta,sum(n.ias)),rep(1:n.cat,n.ias))
	else{
		if(!is.list(beta))
			stop("beta must be either a numeric value or a list.")
		if(any(sapply(beta,length)!=n.ias))
			stop("For each of the interactions in list.ia, a value of beta must be specified.")
	}
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
	getIA<-function(snp,ia) 
		paste("(SNP",snp,ifelse(ia<0," != "," == "),abs(ia),")",sep="",collapse=" & ")
	vec.models<-vec.ias<-character(n.cat)
	for(i in 1:n.cat){
		if(n.ias[i]==1)
			tmp<-getIA(list.snp[[i]],list.ia[[i]])
		else{
			tmp<-character(n.ias[i])
			for(j in 1:n.ias[i])
				tmp[j]<-getIA(list.snp[[i]][[j]],list.ia[[i]][[j]])
		}
		vec.ias[i]<-paste("(",tmp,")",sep="",collapse=" | ")
		vec.models[i]<-paste(beta0[i],"+",paste(beta[[i]]," * (",tmp,")",sep="",
			collapse=" + "))
	}
	mat<-buildSNPmat(n.obs,n.snp,maf,vec.ias)
	mat.prob <- matrix(0,n.obs,n.cat)
	mat2<-as.data.frame(mat)
	for(i in 1:n.cat)
		mat.prob[,i]<-with(mat2, eval(parse(text=vec.models[i])))
	ref.out <- if(withRef) getResponseRef(mat.prob, n.cat, n.obs, beta0, sample.y=sample.y)
		else getResponseCat1(mat.prob, n.cat, n.obs, beta0, n.ias, sample.y=sample.y)
	tab <- table(ref.out$vec.which, ref.out$cl)
	tab <- matrix(tab, ncol=n.cat+withRef, dimnames=dimnames(tab))
	tab.explain<-data.frame(tab,IA=c("None",vec.ias),check.names=FALSE,stringsAsFactors=FALSE)
	out<-list(x=mat,y=ref.out$cl,models=vec.models,maf=maf,tab.explain=tab.explain)
	class(out)<-"simSNPcatResponse"
	out		
}

