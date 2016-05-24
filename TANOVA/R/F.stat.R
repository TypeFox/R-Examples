F.stat<-function(data,f1,f2,type,equal.size=FALSE,trim=0,eb=FALSE){
	if (length(f2)==1){
		z<-matrix(nrow=dim(data)[1],ncol=dim(data)[2])
		z12<-matrix(nrow=dim(data)[1],ncol=dim(data)[2])
		v<-rep(apply(data,1,mean),times=dim(data)[2])
		z<-matrix(v,nrow=dim(data)[1],byrow=FALSE)
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
		S1<-z12-z
		S2<-data-z12
		F<-vector()
		P<-S1^2
		Q<-S2^2
		if (eb==TRUE){
			pr<-prior.sigma(Q,f1,f2)
			F<-apply(P,1,mean,trim=trim)/((apply(Q,1,mean,trim=trim)*length(f1)+pr$v0*pr$s0)/(pr$df+pr$v0))
		}
		if (eb==FALSE){
			F<-apply(P,1,mean,trim=trim)/apply(Q,1,mean,trim=trim)
		}
		return(F)
	}
	if(length(f2)>1){
		g.ix<-group.ix(f1,f2)
		z<-matrix(nrow=dim(data)[1],ncol=dim(data)[2])
		z12<-matrix(nrow=dim(data)[1],ncol=dim(data)[2])
		n1<-nlevels(as.factor(f1))
		n2<-nlevels(as.factor(f2))
		for (i in 1:(n1*n2)){
			ix<-g.ix[[i]]
			if (length(ix)>1){
				v<-rep(apply(data[,ix],1,mean,trim=trim),times=length(ix))
			}
			if (length(ix)==1){
				v<-data[,ix]
			}
			z12[,ix]<-v
		}
		v<-rep(apply(data,1,mean,trim=trim),times=dim(data)[2])
		z<-matrix(v,nrow=dim(data)[1],byrow=FALSE)
		
		if (type==2){
			S1<-z12-z
			S2<-data-z12
		}
		
		if (type!=2){
			if(equal.size==TRUE){
				z1<-matrix(nrow=dim(data)[1],ncol=dim(data)[2])
				z2<-matrix(nrow=dim(data)[1],ncol=dim(data)[2])
				for (i in 1:n1){
					ix<-which(f1==i)
					if (length(ix)>1){
						v<-rep(apply(data[,ix],1,mean,trim=trim),times=length(ix))
					}
					if (length(ix)==1){
						v<-data[,ix]
					}
					z1[,ix]<-v
				}
				for (i in 1:n2){
					ix<-which(f2==i)
					if (length(ix)>1){
						v<-rep(apply(data[,ix],1,mean,trim=trim),times=length(ix))
					}
					if (length(ix)==1){
						v<-data[,ix]
					}
					z2[,ix]<-v
				}
				if (type==1){
					S1<-z12-z1-z2+z
				}
				if (type==3){
					S1<-z1-z
				}
				if (type==4){
					S1<-z2-z
				}
				S2<-data-z12
			}
			
			if (equal.size==FALSE){
				y12<-matrix(nrow=dim(data)[1],ncol=n1*n2)
				y1<-matrix(nrow=dim(data)[1],ncol=n1*n2)
				y2<-matrix(nrow=dim(data)[1],ncol=n1*n2)
				y<-matrix(nrow=dim(data)[1],ncol=n1*n2)
				for (i in 1:(n1*n2)){
					ix<-g.ix[[i]]
					if (length(ix)>1){
						v<-apply(data[,ix],1,mean,trim=trim)
					}
					if (length(ix)==1){
						v<-data[,ix]
					}
					y12[,i]<-v
				}
				g1<-rep(c(1:n1),each=n2)
				g2<-rep(c(1:n2),times=n1)
				for (i in 1:n1){
					ix<-which(g1==i)
					v<-rep(apply(y12[,ix],1,mean,trim=0),times=length(ix))
					y1[,ix]<-v
				}
				for (i in 1:n2){
					ix<-which(g2==i)
					v<-rep(apply(y12[,ix],1,mean,trim=0),times=length(ix))
					y2[,ix]<-v
				}
				v<-rep(apply(y12,1,mean,trim=0),each=n1*n2)
				y<-matrix(v,nrow=dim(data)[1],byrow=TRUE)
				if (type==1){
					S1<-y12-y1-y2+y
				}
				if (type==3){
					S1<-y1-y
				}
				if (type==4){
					S1<-y2-y
				}
				S2<-data-z12
			}
		}
		F<-vector()
		P<-S1^2
		Q<-S2^2
		if (eb==TRUE){
			pr<-prior.sigma(Q,f1,f2)
		}
		if (equal.size==TRUE){
			if (eb==TRUE){
				F<-apply(P,1,mean,trim=trim)/((apply(Q,1,mean,trim=trim)*length(f1)+pr$v0*pr$s0)/(pr$df+pr$v0))
			}
			if (eb==FALSE){
				F<-apply(P,1,mean,trim=trim)/apply(Q,1,mean,trim=trim)
			}
		}
		if (equal.size==FALSE){
			if (eb==TRUE){
				F<-apply(P,1,mean)/((apply(Q,1,mean,trim=trim)*length(f1)+pr$v0*pr$s0)/(pr$df+pr$v0))
			}
			if (eb==FALSE){
				F<-apply(P,1,mean)/apply(Q,1,mean,trim=trim)
			}
		}
		return(F)
	}
}