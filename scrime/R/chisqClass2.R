`chisqClass2` <-
function(data,cl,n.cat,same.ncat=FALSE,compPval=TRUE){
	if(any(is.na(cl)))
		stop("No missing values allowed in cl.")
	if(missing(n.cat))
		n.cat<-max(data,na.rm=TRUE)
	uni.cl<-sort(unique(cl))
	n.lev<-length(uni.cl)
	if(length(n.lev)>10)
		stop("cl contains more than 10 different values.")
	if(any(uni.cl!=1:n.lev))
		stop("The labels of the classes must be 1 to ",n.lev,".")
	n.obs<-ncol(data)
	if(length(cl)!=n.obs)
		stop("The length of cl must be equal to the number of observations.")
	CL<-matrix(0,n.obs,n.lev)
	for(i in 1:n.lev)
		CL[cl==i,i]<-1
	if(any(colSums(CL)<2))
		stop("There must be at least two observations per class.")
	listCells<-vector("list",n.cat)
	stats<-numeric(nrow(data))
	if(any(is.na(data))){
		anyna<-TRUE
		mat.na<-!is.na(data)
		vec.ncl<-mat.na%*%CL
		n.obs<-rowSums(mat.na)
	}
	else{
		mat.na<-TRUE
		vec.ncl<-colSums(CL)
		anyna<-FALSE
	}
	if(compPval)
		df<-numeric(nrow(data))
	for(i in 1:n.cat){
		mat.i<-mat.na & data==i
		N<-mat.i%*%CL
		rS<-rowSums(N)
		if(same.ncat && any(rS==0))
			stop("All rows of data must show the same number of categories.",call.=FALSE)
		tmp<-if(anyna) rS*vec.ncl/n.obs  else rS%*%t(vec.ncl)/n.obs
		tmp[tmp==0]<-1
		tmp2<-rowSums(N*N/tmp)
		stats<-tmp2+stats
		if(compPval)
			df <- df + (rowSums(mat.i)>0)
	}
	names(stats)<-rownames(data)
	stats<-stats-n.obs
	if(!compPval)
		return(stats)
	df<-(df-1)*(n.lev-1)
	rawp<-pchisq(stats,df,lower.tail=FALSE)
	rawp[df==0] <- 1
	structure(list(stats=stats,df=df,rawp=rawp))
}

