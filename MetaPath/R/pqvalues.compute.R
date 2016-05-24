pqvalues.compute <-
function(Stat.0,Stat.B,Stat.type,PI0=1){


## Stat.0 GX1 matrix
## Stat.B GXB matrix

	Stat.0=as.matrix(abs(Stat.0))
	Stat.B=as.matrix(abs(Stat.B))
	
	if (nrow(Stat.0)!=nrow(Stat.B)) stop('# of rows of Stat.0 and Stat.B should be same')

	B=ncol(Stat.B)
	G=nrow(Stat.B)

	if (Stat.type=='Tstat'){

	Stat.all=cbind(Stat.0,Stat.B)
	count.0=apply(Stat.0,1,function(x) sum(x<=Stat.0,na.rm=T))
	Stat.rank.all=rank(-as.vector(Stat.all),ties.method= "max")
	pvalue.0=as.matrix(Stat.rank.all[1:G]-count.0)/(B*G)

	if (is.null(PI0)){
	PI0=sum(pvalue.0>=.5)/(.5*G)
	} else {
	PI0=1 }
	
	
	qvalue.0=PI0*pvalue.0*G/(apply(Stat.0,1,function(x) sum(x<=Stat.0,na.rm=T)))
	qvalue.0=ifelse(qvalue.0<=1,qvalue.0,1)

	Stat.rank=rank(-as.vector(Stat.B),ties.method= "max")
	pvalue.B=matrix(Stat.rank/(B*G),G,B)
	rownames(pvalue.B)=rownames(Stat.B)
	colnames(pvalue.B)=colnames(Stat.B)

	} else if (Stat.type=='Pvalue'){

	Stat.all=cbind(Stat.0,Stat.B)
	count.0=apply(Stat.0,1,function(x) sum(x>=Stat.0))
	Stat.rank.all=rank(as.vector(Stat.all),ties.method= "max")
	pvalue.0=as.matrix(Stat.rank.all[1:G]-count.0)/(B*G)

	if (is.null(PI0)){
	PI0=sum(pvalue.0>=.5)/(.5*G)
	} else {
	PI0=1 }
	
	
	qvalue.0=PI0*pvalue.0*G/(apply(Stat.0,1,function(x) sum(x>=Stat.0,na.rm=T)))
	qvalue.0=ifelse(qvalue.0<=1,qvalue.0,1)

	Stat.rank=rank(as.vector(Stat.B),ties.method= "max")
	pvalue.B=matrix(Stat.rank/(B*G),G,B)
	rownames(pvalue.B)=rownames(Stat.B)
	colnames(pvalue.B)=colnames(Stat.B)

	} else {
		stop('wrong stat.type')
	}
	
	return(list(pvalue.0=pvalue.0,pvalue.B=pvalue.B,qvalue.0=qvalue.0,PI0=PI0))
		
}
