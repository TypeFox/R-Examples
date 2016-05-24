F.stat2<-function(data,f1,f2,tp,type,trim=0,eb=FALSE){
	
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
