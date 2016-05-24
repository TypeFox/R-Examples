`pamCat` <-
function(data,cl,theta=NULL,n.theta=10,newdata=NULL,newcl=NULL){
	if(is.null(rownames(data)))
		stop("data must have row names.")
	n.cat<-max(data)
	out<-computeContClass(data,cl,n.cat)
	N<-out$mat.obs
	Nexp<-out$mat.exp
	rm(out)
	n.lev<-length(unique(cl))
	vec.ncl<-table(cl)
	mat.chisq<-matrix(0,nrow(data),n.lev)
	for(i in 1:n.lev){
		ids<-i+n.lev*(0:(n.cat-1))
		mat.chisq[,i]<-rowSums(N[,ids]*N[,ids]/Nexp[,ids])-vec.ncl[i]
	}
	rownames(mat.chisq)<-rownames(data)
	colnames(mat.chisq)<-names(vec.ncl)
	if(is.null(theta)){
		if(n.theta<1)
			stop("n.theta must be at least 1.")
		rangeChisq<-range(mat.chisq)
		theta<-seq(max(0.1,rangeChisq[1]),max(1,rangeChisq[2]-0.1),le=n.theta)
		theta<-unique(round(theta,1))
		n.theta<-length(theta)
	}
	else{
		if(any(theta<=0))
			stop("theta must be larger than 0.")
		n.theta<-length(theta)
	}
	if(is.null(newdata))
		newdata<-data
	else{
		if(is.null(newcl))
			stop("If newdata is specified, newcl must also be specified.")
	}
	if(is.null(newcl))
		newcl<-cl
	if(ncol(newdata)!=length(newcl))
		stop("The length of newcl must be equal to the number of columns of newdata.")
	mat.theta<-cbind(Theta=theta,Variables=0,MCR=0)
	out<-list(mat.chisq=mat.chisq,mat.obs=N,mat.exp=Nexp,mat.theta=mat.theta,tab.cl=vec.ncl,n.cat=n.cat)
	class(out)<-"pamCat"
	for(i in 1:n.theta){
		tmp<-predict(out,newdata,theta=theta[i],add.nvar=TRUE)
		mat.theta[i,2]<-tmp$n.var
		mat.theta[i,3]<-mean(tmp$pred!=newcl)
	}
	out$mat.theta<-as.data.frame(mat.theta)
	out	
}

