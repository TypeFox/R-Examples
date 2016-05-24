`computeContClass` <-
function(data,cl,n.cat){
	if(!is.matrix(data))
		stop("data must be a matrix.",call.=FALSE)
	if(any(is.na(cl)))
		stop("No missing values allowed in cl.",call.=FALSE)
	if(any(is.na(data)))
		stop("No missing values allowed in data.",call.=FALSE)
	n.obs<-ncol(data)
	if(length(cl)!=n.obs)
		stop("The length of cl must be equal to the number of columns of data.",call.=FALSE)
	uni.cl<-sort(unique(cl))
	n.lev<-length(uni.cl)
	if(length(n.lev)>10)
		stop("cl contains more than 10 different values.",call.=FALSE)
	if(any(uni.cl!=1:n.lev))
		stop("The labels of the classes must be 1 to ",n.lev,".",call.=FALSE)
	if(missing(n.cat))
		n.cat<-max(data)
	checkCatMat(data,n.cat)
	CL<-matrix(0,n.obs,n.lev)
	for(i in 1:n.lev)
		CL[cl==i,i]<-1
	vec.ncl<-colSums(CL)
	if(any(vec.ncl<2))
		stop("There must be at least two observations per class.")
	mat.obs<-mat.exp<-matrix(0,nrow(data),n.cat*n.lev)
	for(i in 1:n.cat){
		tmp<-data==i
		rs<-rowSums(tmp)
		if(any(rs==0))
			stop("All variables must show the same number of categories.",call.=FALSE)
		ids<-(i-1)*n.lev+(1:n.lev)
		mat.obs[,ids]<-tmp%*%CL
		mat.exp[,ids]<-rs%*%t(vec.ncl)/n.obs
	}
	rownames(mat.obs)<-rownames(mat.exp)<-rownames(data)
	colnames(mat.obs)<-colnames(mat.exp)<-paste("N",rep(1:n.lev,n.cat),rep(1:n.cat,e=n.lev),sep="")
	return(list(mat.obs=mat.obs,mat.exp=mat.exp))
}

