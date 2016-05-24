NANOVA.test2<-function(data,f1,f2,type,time.course,equal.size=FALSE,B=100,robustify=FALSE,eb=FALSE,df=0){
	if (robustify==FALSE){
		tm<-0
	}
	if (robustify==TRUE){
		tm<-0.2
	}
	a<-proj.dir(data=data,f1=f1,f2=f2,time.course=time.course,type=type,df=df)
	d<-proj.data(data=data,time.course=time.course,a=a)
	F<-F.stat(data=d,f1=f1,f2=f2,type=type,equal.size=equal.size,trim=tm,eb=eb)
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
		F.null<-matrix(nrow=dim(data)[1]/time.course,ncol=B)
		for (i in 1:B){
			ix<-sample(c(1:p),replace=TRUE,size=p)
			d.null<-M+r[,ix]
			a<-proj.dir(data=d.null,f1=f1,f2=f2,time.course=time.course,type=type,df=df)
			d<-proj.data(data=d.null,time.course=time.course,a=a)
			F.null[,i]<-F.stat(data=d,f1=f1,f2=f2,type=type,trim=tm,equal.size=equal.size,eb=eb)
		}
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
		F.null<-matrix(nrow=dim(data)[1]/time.course,ncol=B)
		for (i in 1:B){
			ix<-sample(c(1:p),replace=TRUE,size=p)
			d.null<-M+r[,ix]
			a<-proj.dir(data=d.null,f1=f1,f2=f2,time.course=time.course,type=type,df=df)
			d<-proj.data(data=d.null,time.course=time.course,a=a)
			F.null[,i]<-F.stat(data=d,f1=f1,f2=f2,type=type,trim=tm,equal.size=equal.size,eb=eb)
		}
	}
	delta<-z.score(F,F.null)$z
	gene<-sort(delta,decreasing=TRUE,index=TRUE)$ix
	p<-vector()
	n<-length(F)
	for (i in 1:n){
		p[i]<-1-ecdf(F.null[i,])(F[i])
}
list(gene.order=gene, F=F, F.null=F.null, pvalue=p, delta=delta)
}