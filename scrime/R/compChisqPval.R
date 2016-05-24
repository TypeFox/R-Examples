`compChisqPval` <-
function(data,stats,n.cat,asMatrix=TRUE){
	n.row<-nrow(data)
	df<-numeric(n.row)
	for(i in 1:n.cat){
		tmp<-rowSums(data==i,na.rm=TRUE)>0
		df<-df+tmp
	}
	df<-df-1
	mat.df<-df%*%t(df)
	df<-mat.df[lower.tri(mat.df)]
	rawp<-pchisq(stats,df,lower.tail=FALSE)
	rawp[df==0] <- 1
	if(!asMatrix)
		return(list(stats=stats,df=df,rawp=rawp))
	mat.stat<-mat.df<-mat.p<-matrix(0,n.row,n.row)
	mat.stat[lower.tri(mat.stat)]<-stats
	mat.df[lower.tri(mat.df)]<-df
	mat.p[lower.tri(mat.p)]<-rawp
	colnames(mat.stat)<-rownames(mat.stat)<-rownames(data)
	colnames(mat.df)<-rownames(mat.df)<-rownames(data)
	colnames(mat.p)<-rownames(mat.p)<-rownames(data)
	return(list(stats=mat.stat,df=mat.df,rawp=mat.p))
}

