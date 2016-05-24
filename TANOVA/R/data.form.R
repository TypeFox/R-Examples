data.form<-function(data,f1,f2,tp){
	if (length(f2)>1){
		n1<-nlevels(as.factor(f1))
		n2<-nlevels(as.factor(f2))
		tps<-length(unique(tp))
		n<-dim(data)[1]
		id<-list()
		l<-1
		k<-1
		fc1<-vector()
		fc2<-vector()
		for (i in 1:n1){
			for (j in 1:n2){
				id[[l]]<-list()
				ix<-which(f1==i&f2==j&tp==1)
				id[[l]]<-seq(from=k,by=1,length.out=length(ix))
				k<-k+length(ix)
				l<-l+1
				fc1<-append(fc1,rep(i,times=length(ix)))
				fc2<-append(fc2,rep(j,times=length(ix)))
			}
		}
		
		d<-matrix(nrow=n*tps,ncol=length(fc1))
		l<-1
		for (i in 1:n1){
			for (j in 1:n2){
				for (k in 1:tps){
					rix<-seq(from=k, by=tps, length.out=n)
					ix<-which(f1==i &f2==j &tp==k)
					d[rix,id[[l]]]<-data[,ix]
				}
				l<-l+1
			}
		}
		return(list(d=d,fc1=fc1,fc2=fc2))
	}
	if (length(f2)==1){
		n1<-nlevels(as.factor(f1))
		tps<-length(unique(tp))
		n<-dim(data)[1]
		id<-list()
		l<-1
		k<-1
		fc1<-vector()
		for (i in 1:n1){
			id[[l]]<-list()
			ix<-which(f1==i&tp==1)
			id[[l]]<-seq(from=k,by=1,length.out=length(ix))
			k<-k+length(ix)
			l<-l+1
			fc1<-append(fc1,rep(i,times=length(ix)))
		}
		d<-matrix(nrow=n*tps,ncol=length(fc1))
		l<-1
		for (i in 1:n1){
			for (k in 1:tps){
				rix<-seq(from=k, by=tps, length.out=n)
				ix<-which(f1==i&tp==k)
				d[rix,id[[l]]]<-data[,ix]
}
l<-l+1
}
return(list(d=d,fc1=fc1,fc2=0))
}
}