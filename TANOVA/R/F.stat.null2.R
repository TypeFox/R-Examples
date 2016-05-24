F.stat.null2<-function(data,f1,f2,tp,type,B=100,trim=trim,eb=FALSE){
	M<-matrix(nrow=dim(data)[1],ncol=dim(data)[2])
	if (length(f2)==1){
		n1<-nlevels(as.factor(f1))
		tps<-length(unique(tp))
		n<-dim(data)[1]
		y<-matrix(nrow=n,ncol=dim(data)[2])
		for (k in 1:tps){
			ix<-which(tp==k)
			if (length(ix)>1){
				v<-apply(data[,ix],1,mean)
			}
			if (length(ix)==1){
				v<-data[,ix]
			}
			y[,ix]<-rep(v,times=length(ix))
		}
		M<-y
		temp<-matrix(nrow=n,ncol=dim(data)[2])
		for (i in 1:n1){
			for (k in 1:tps){
				ix<-which(f1==i & tp==k)
				if (length(ix)>1){
					v<-apply(data[,ix],1,mean)
				}
				if (length(ix)==1){
					v<-data[,ix]
				}
				temp[,ix]<-rep(v,times=length(ix))
			}
		}
		r<-data-temp
	}
	
	if (length(f2)>1){
		n1<-length(unique(f1))
		n2<-length(unique(f2))
		tps<-length(unique(tp))
		n<-dim(data)[1]
		temp<-matrix(nrow=n,ncol=dim(data)[2])
		for (i in 1:n1){
			for (j in 1:n2){
				for (k in 1:tps){
					ix<-which(f1==i & f2==j & tp==k)
					if (length(ix)>1){
						v<-apply(data[,ix],1,mean)
					}
					if (length(ix)==1){
						v<-data[,ix]
					}
					temp[,ix]<-rep(v,times=length(ix))
				}
			}
		}
		y12<-temp
		x12<-matrix(nrow=n*tps,ncol=n1*n2)
		r1<-rep(c(1:n1),each=n2)
		r2<-rep(c(1:n2),times=n1)
		for (i in 1:n1){
			for (j in 1:n2){
				for (k in 1:tps){
					rix<-seq(from=k,by=tps,length.out=n)
					cix<-which(r1==i & r2==j)
					ix<-which(f1==i & f2==j & tp==k)
					if (length(ix)>1){
						v<-apply(data[,ix],1,mean)
					}
					if (length(ix)==1){
						v<-data[,ix]
					}
					x12[rix,cix]<-v
				}
			}
		}
		
		y<-matrix(nrow=n,ncol=dim(data)[2])
		for (k in 1:tps){
			rix<-seq(from=k, by=tps,length.out=n)
			ix<-which(tp==k)
			y[,ix]<-rep(apply(x12[rix,],1,mean),times=length(ix))
		}
		
		if (type==1){
			y1<-matrix(nrow=n,ncol=dim(data)[2])
			for (i in 1:n1){
				for (k in 1:tps){
					rix<-seq(from=k, by=tps,length.out=n)
					ix<-which(f1==i & tp==k)
					y1[,ix]<-rep(apply(x12[rix,r1==i],1,mean),times=length(ix))
				}
			}
			
			y2<-matrix(nrow=n,ncol=dim(data)[2])
			for (i in 1:n2){
				for (k in 1:tps){
					rix<-seq(from=k, by=tps,length.out=n)
					ix<-which(f2==i & tp==k)
					y2[,ix]<-rep(apply(x12[rix,r2==i],1,mean),times=length(ix))
				}
			}
			M<-y1+y2-y
		}
		
		if (type==2){
			M<-y
		}
		
		if(type==3){
			y2<-matrix(nrow=n,ncol=dim(data)[2])
			for (i in 1:n2){
				for (k in 1:tps){
					rix<-seq(from=k, by=tps,length.out=n)
					ix<-which(f2==i & tp==k)
					y2[,ix]<-rep(apply(x12[rix,r2==i],1,mean),times=length(ix))
				}
			}
			M<-y2
		}
		
		if (type==4){
			y1<-matrix(nrow=n,ncol=dim(data)[2])
			for (i in 1:n1){
				for (k in 1:tps){
					rix<-seq(from=k, by=tps,length.out=n)
					ix<-which(f1==i & tp==k)
					y1[,ix]<-rep(apply(x12[rix,r1==i],1,mean),times=length(ix))
				}
			}
			M<-y1
		}
		r<-data-temp
	}
	F.null<-matrix(nrow=dim(data)[1],ncol=B)
	p<-dim(data)[2]
	for (i in 1:B){
		ix<-sample(c(1:p),replace=TRUE,size=p)
		d.null<-M+r[,ix]
		F.null[,i]<-F.stat2(d.null,f1,f2,tp=tp,type=type,trim=trim,eb=eb)
	}
	return(F.null)
}