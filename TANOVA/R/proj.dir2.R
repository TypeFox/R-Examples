proj.dir2<-function(data,f1,f2,tp,type=2,trim=0,...){
	if (length(f2)==1){
		n1<-nlevels(as.factor(f1))
		tps<-nlevels(as.factor(tp[f1==1]))
		n<-dim(data)[1]
		a<-matrix(nrow=n,ncol=tps)
		u<-matrix(nrow=n*n1,ncol=tps)
		for (i in 1:n1){
			z<-(c(1:n)-1)*n1+i
			for (j in 1:tps){
				ix<-which(f1==i & tp==j)
				if (length(ix)>1){
					v<-apply(data[,ix],1,mean,trim=trim)
				}
				if (length(ix)==1){
					v<-data[,ix]
				}
				u[z,j]<-v
			}
		}
		e<-matrix(1,nrow=n1,ncol=1)
		for (i in 1:n){
			z<-c(((i-1)*n1+1):(i*n1))
			a[i,]<-Re(eigen(t(u[z,])%*%u[z,]-1/n1*t(u[z,])%*%e%*%t(e)%*%u[z,])$vectors[,1])
			a[i,]<-sign(a[i,1])*a[i,]
		}
		return(a)
	}
	
	if (length(f2)>1){
		n1<-nlevels(as.factor(f1))
		n2<-nlevels(as.factor(f2))
		tps<-length(unique(tp))
		n<-dim(data)[1]
		y12<-matrix(nrow=n1*n2*n,ncol=tps)
		y<-matrix(nrow=n1*n2*n,ncol=tps)
		y1<-matrix(nrow=n1*n2*n,ncol=tps)
		y2<-matrix(nrow=n1*n2*n,ncol=tps)
		r1<-rep(rep(c(1:n1),each=n2),times=n)
		r2<-rep(rep(c(1:n2),times=n1),times=n)
		for (i in 1:n1){
			for (j in 1:n2){
				for (k in 1:tps){
					ix<-which(f1==i &  f2==j & tp==k)
					if (length(ix)>1){
						v<-apply(data[,ix],1,mean,trim=trim)
					}
					if (length(ix)==1){
						v<-data[,ix]
					}
					y12[which(r1==i & r2==j),k]<-v
				}
			}
		}
		
		for (i in 1:n1){
			for (k in 1:tps){
				temp<-matrix(y12[r1==i,k],nrow=n,byrow=TRUE)
				y1[r1==i,k]<-rep(apply(temp,1,mean,trim=trim),each=n2)
			}
		}
		
		for (j in 1:n2){
			for (k in 1:tps){
				temp<-matrix(y12[r2==j,k],nrow=n,byrow=TRUE)
				y2[r2==j,k]<-rep(apply(temp,1,mean,trim=trim),each=n1)
			}
		}
		
		
		for (k in 1:tps){
			temp<-matrix(y12[,k],nrow=n,byrow=TRUE)
			y[,k]<-rep(apply(temp,1,mean,trim=trim),each=n1*n2)
		}
		
		if (type==1){
			X<-y12-y1-y2+y
		}
		if (type==2){
			X<-y12-y
		}
		if (type==3){
			X<-y1-y
		}
		if (type==4){
			X<-y2-y
		}
		a<-matrix(nrow=n,ncol=tps)
		for (i in 1:n){
			ix<-c(((i-1)*n1*n2+1):(i*n1*n2))
			x<-X[ix,]
			a[i,]<-Re(eigen(t(x)%*%(x))$vector[,1])
			a[i,]<-sign(a[i,1])*a[i,]
		}
		return(a)
	}
}

F.stat2<-function(data,f1,f2,tp,type,trim=trim,eb=FALSE){
	if (length(f2)==1){
		n1<-nlevels(as.factor(f1))
		tps<-nlevels(as.factor(tp[f1==1]))
		n<-dim(data)[1]
		u<-matrix(nrow=tps*n,ncol=n1)
		for (i in 1:n1){
			for (j in 1:tps){
				z<-seq(from=j,by=tps,length.out=n)
				ix<-which(f1==i & tp==j)
				if (length(ix)>1){
					v<-apply(data[,ix],1,mean,trim=trim)
				}
				if (length(ix)==1){
					v<-data[,ix]
				}
				u[z,i]<-v
			}
		}
		y<-matrix(rep(apply(u,1,mean,trim=trim),each=n1),nrow=n*tps,byrow=TRUE)
		S1<-u-y
		temp<-matrix(nrow=n,ncol=dim(data)[2])
		for (i in 1:n1){
			for (k in 1:tps){
				ix<-which(f1==i & tp==k)
				if (length(ix)>1){
					v<-apply(data[,ix],1,mean,trim=trim)
				}
				if (length(ix)==1){
					v<-data[,ix]
				}
				temp[,ix]<-matrix(rep(v,each=length(ix)),nrow=n,byrow=TRUE)
			}
		}
		S2<-data-temp
		if (eb==TRUE){
			pr<-prior.sigma(S2^2,f1,f2,tp)
		}
		F<-vector()
		a<-proj.dir2(data,f1,f2,tp,trim=trim)
		for (i in 1:n){
			ix<-c(((i-1)*tps+1):(i*tps))
			B<-S1[ix,]%*%t(S1[ix,])
			if (eb==TRUE){
				F[i]<-mean(a[i,]%*%B*a[i,],trim=trim)/((mean(S2[i,]^2,trim=trim)*length(f1)+pr$v0*pr$s0)/(pr$df+pr$v0))
			}
			if (eb==FALSE){
				F[i]<-mean(a[i,]%*%B*a[i,],trim=trim)/mean(S2[i,]^2,trim=trim)
			}
		}
		return(F)
	}
	if (length(f2)>1){
		n1<-length(unique(f1))
		n2<-length(unique(f2))
		tps<-length(unique(tp))
		n<-dim(data)[1]
		
		y12<-matrix(nrow=n*tps,ncol=n1*n2)
		r1<-rep(c(1:n1),each=n2)
		r2<-rep(c(1:n2),times=n1)
		for (i in 1:n1){
			for (j in 1:n2){
				for (k in 1:tps){
					rix<-seq(from=k,by=tps,length.out=n)
					cix<-which(r1==i & r2==j)
					ix<-which(f1==i & f2==j & tp==k)
					if (length(ix)>1){
						v<-apply(data[,ix],1,mean,trim=trim)
					}
					if (length(ix)==1){
						v<-data[,ix]
					}
					y12[rix,cix]<-v
				}
			}
		}
		v<-rep(apply(y12,1,mean,trim=trim),each=n1*n2)
		y<-matrix(v,nrow=n*tps,byrow=TRUE)
		
		if (type==1){
			y1<-matrix(nrow=n*tps,ncol=n1*n2)
			for (i in 1:n1){
				ix<-which(r1==i)
				y1[,ix]<-rep(apply(y12[,ix],1,mean,trim=trim),times=length(ix))
			}
			
			y2<-matrix(nrow=n*tps,ncol=n1*n2)
			for (j in 1:n2){
				ix<-which(r2==j)
				y2[,ix]<-rep(apply(y12[,ix],1,mean,trim=trim),times=length(ix))
			}
			
			S1<-y12-y1-y2+y
		}
		
		if (type==2){
			S1<-y12-y
		}
		
		if(type==3){
			y1<-matrix(nrow=n*tps,ncol=n1*n2)
			for (i in 1:n1){
				ix<-which(r1==i)
				y1[,ix]<-rep(apply(y12[,ix],1,mean,trim=trim),times=length(ix))
			}
			S1<-y1-y
		}
		
		if (type==4){
			y2<-matrix(nrow=n*tps,ncol=n1*n2)
			for (j in 1:n2){
				ix<-which(r2==j)
				y2[,ix]<-rep(apply(y12[,ix],1,mean,trim=trim),times=length(ix))
			}
			S1<-y2-y
		}
		temp<-matrix(nrow=n,ncol=dim(data)[2])
		for (i in 1:n1){
			for (j in 1:n2){
				for (k in 1:tps){
					ix<-which(f1==i & f2==j & tp==k)
					if (length(ix)>1){
						v<-apply(data[,ix],1,mean,trim=trim)
					}
					if (length(ix)==1){
						v<-data[,ix]
					}
					temp[,ix]<-rep(v,times=length(ix))
				}
			}
		}
		
		S2<-data-temp
#sigma2<-apply(S2^2,1,mean)
		if (eb==TRUE){
			pr<-prior.sigma(S2^2,f1,f2,tp)
		}
		F<-vector()
		a<-proj.dir2(data,f1,f2,tp,type=type,trim=trim)
		for (i in 1:n){
			ix<-c(((i-1)*tps+1):(i*tps))
			B<-S1[ix,]%*%t(S1[ix,])
			if (eb==TRUE){
				F[i]<-mean(a[i,]%*%B*a[i,],trim=trim)/((mean(S2[i,]^2,trim=trim)*length(f1)+pr$v0*pr$s0)/(pr$df+pr$v0))
			}
			if (eb==FALSE){
				F[i]<-mean(a[i,]%*%B*a[i,],trim=trim)/mean(S2[i,]^2,trim=trim)
			}
		}
		return(F)
	}
}


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
