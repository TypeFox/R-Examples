F.stat.null<-function(data,f1,f2,type,trim=0,B=100,equal.size=FALSE,eb=FALSE){
	if (length(f2)==1){
		M<-matrix(nrow=dim(data)[1],ncol=dim(data)[2])
		M<-rep(apply(data,1,mean), times=length(f1))
		z12<-matrix(nrow=dim(data)[1],ncol=dim(data)[2])
		n1<-nlevels(as.factor(f1))
		for (i in 1:n1){
			ix<-which(f1==i)
			if (length(ix)>1){
				v<-rep(apply(data[,ix],1,mean),times=length(ix))
			}
			if (length(ix)==1){
				v<-data[,ix]
			}
			z12[,ix]<-v
		}
		r<-data-z12
		p<-length(f1)
		F.null<-matrix(nrow=dim(data)[1],ncol=B)
		for (i in 1:B){
			ix<-sample(c(1:p),replace=TRUE,size=p)
			d.null<-M+r[,ix]
			F.null[,i]<-F.stat(d.null,f1,f2,type=type,trim=trim,equal.size=equal.size,eb=eb)
		}
		return(F.null)
	}
	if (length(f2)>1){
		if (type==1){
			m<-ls.estimate(data,f1,f2)
			M<-m$Mab
		}
		if (type==2){
			m<-ls.estimate(data,f1,f2)
			M<-m$M0
		}
		if (type==3){
			m<-ls.estimate(data,f1,f2)
			M<-m$Mb
		}
		if (type==4){
			m<-ls.estimate(data,f2,f1)
			M<-m$Mb
		}
		p<-length(f1)
		g.ix<-group.ix(f1,f2)
		
		z12<-matrix(nrow=dim(data)[1],ncol=dim(data)[2])
		n1<-nlevels(as.factor(f1))
		n2<-nlevels(as.factor(f2))
		for (i in 1:(n1*n2)){
			ix<-g.ix[[i]]
			if (length(ix)>1){
				v<-rep(apply(data[,ix],1,mean),times=length(ix))
			}
			if (length(ix)==1){
				v<-data[,ix]
			}
			z12[,ix]<-v
		}
		r<-data-z12
		F.null<-matrix(nrow=dim(data)[1],ncol=B)
		for (i in 1:B){
			ix<-sample(c(1:p),replace=TRUE,size=p)
			d.null<-M+r[,ix]
			F.null[,i]<-F.stat(d.null,f1,f2,type=type,trim=trim,equal.size=equal.size,eb=eb)
}
return(F.null)
}
}