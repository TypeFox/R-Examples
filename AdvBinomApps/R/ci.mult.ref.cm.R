ci.mult.ref.cm <- function(k,n,A.ref,A.follow,K,theta,alpha=0.1,p.target=1,prec=2,tailcut=1e-08,tol=1e-12){
	
	if((alpha<=0)||(alpha>=1)||(p.target<=0)||(p.target>1)){
		return<-"Input error: check alpha and/or p.target!"
	}else if(A.follow<=0){
		return<-"Input error: A.follow has to be larger than 0!"
	}else{
		scale<-phi.mult.ref.cm(k,n,A.ref,K,theta,prec,tailcut)
		if(is.character(scale)){
			return<-scale
		}else{
			A.prec<-round(A.ref,prec)	
			A.prec[A.prec-A.ref>0]<-A.prec[A.prec-A.ref>0]-10^(-prec)
			A.ref<-A.prec	

			phi<-scale$phi.cm
			A.gcd<-scale$A.gcd
			r<-length(k)

			n.gcd<-sum(n*A.ref/A.gcd)
      		f.p.gcd.cm<-function(p.gcd.cm){as.numeric(phi$prob%*%pbinom(phi$k.gcd,n.gcd,p.gcd.cm))-alpha}
     			p.gcd.cm.hat<-uniroot(f.p.gcd.cm,lower=0,upper=1,tol=tol)$root
      
			p.ref.cm<-NULL  
			for(j in c(1:r)){
				if(k[j]>0){
					p.ref.cm<-c(p.ref.cm,cm.clopper.pearson.ci(n[j],c(K[j,],k[j]-sum(K[j,])),c(theta,0),alpha,uniroot.tol=tol)$Upper.limit)
				}else{
					p.ref.cm<-c(p.ref.cm,qbeta(1-alpha,1,n[j]))
				}
			}

			p.mm.cm<-1-(1-p.gcd.cm.hat)^(1/A.gcd)
			p.follow.cm<-1-(1-p.mm.cm)^A.follow

			if(p.follow.cm>p.target){	
				n.add<-NULL
				for(i in c(1:r)){
					f.n<-function(n.add.i){
						vec<-rep(0,r);vec[i]<-floor(n.add.i)
						n.new<-n+vec
						scale<-phi.mult.ref.cm(k,n.new,A.ref,K,theta,prec,tailcut)
						phi<-scale$phi.cm
						A.gcd<-scale$A.gcd
						n.gcd<-sum(n.new*A.ref/A.gcd)
      					f.p.gcd.cm<-function(p.gcd.cm){as.numeric(phi$prob%*%pbinom(phi$k.gcd,n.gcd,p.gcd.cm))-alpha}
     						p.gcd.cm.hat<-uniroot(f.p.gcd.cm,lower=0,upper=1,tol=tol)$root
						p.fol<-1-(1-p.gcd.cm.hat)^(A.follow/A.gcd)
						p.fol-p.target
					}
					n.add.i<-ceiling(uniroot(f.n,lower=0,upper=100*n[i],extendInt="yes")$root)
					n.add<-c(n.add,n.add.i)
				}
			}else{
				n.add<-0
			}
			return<-list(p.ref.cm=p.ref.cm,p.mm.cm=p.mm.cm,p.follow.cm=p.follow.cm,n.add.cm=n.add)
		}
	}
	return
}
