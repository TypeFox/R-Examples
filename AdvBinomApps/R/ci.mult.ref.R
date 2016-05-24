ci.mult.ref <- function(k,n,A.ref,A.follow,alpha=0.1,p.target=1,prec=2,tailcut=1e-08,tol=1e-12){

	if((alpha<=0)||(alpha>=1)||(p.target<=0)||(p.target>1)){
		return<-"Input error: check alpha and/or p.target!"
	}else if(A.follow<=0){
		return<-"Input error: A.follow has to be larger than 0!"
	}else{

		r<-length(k)
		scale<-phi.mult.ref(k,n,A.ref,prec,tailcut)

		if(is.character(scale)){
			return<-scale
		}else{
			phi<-scale$phi
			A.gcd<-scale$A.gcd

			A.prec<-round(A.ref,prec)	
			A.prec[A.prec-A.ref>0]<-A.prec[A.prec-A.ref>0]-10^(-prec)
			A.ref<-A.prec

			n.gcd<-sum(n*A.ref/A.gcd)
     			f.p.gcd<-function(p.gcd){as.numeric(phi$prob%*%pbinom(phi$k.gcd,n.gcd,p.gcd))-alpha}
  			p.gcd.hat<-uniroot(f.p.gcd,lower=0,upper=1,tol=tol)$root
        
			p.ref<-qbeta(1-alpha,k+1,n-k)
			p.mm<-1-(1-p.gcd.hat)^(1/A.gcd)
			p.follow<-1-(1-p.mm)^A.follow

			if(p.follow>p.target){	
				n.add<-NULL
				for(i in c(1:r)){
					f.n<-function(n.add.i){
						vec<-rep(0,r);vec[i]<-floor(n.add.i)
						n.new<-n+vec
						scale<-phi.mult.ref(k,n.new,A.ref,prec,tailcut)
						phi<-scale$phi
						A.gcd<-scale$A.gcd
						n.gcd<-sum(n.new*A.ref/A.gcd)
      					f.p.gcd<-function(p.gcd){as.numeric(phi$prob%*%pbinom(phi$k.gcd,n.gcd,p.gcd))-alpha}
     						p.gcd.hat<-uniroot(f.p.gcd,lower=0,upper=1,tol=tol)$root
						p.fol<-1-(1-p.gcd.hat)^(A.follow/A.gcd)
						p.fol-p.target
					}
					n.add.i<-ceiling(uniroot(f.n,lower=0,upper=100*n[i],extendInt="yes")$root)
					n.add<-c(n.add,n.add.i)
				}
			}else{
				n.add<-0
			}
			return<-list(p.ref=p.ref,p.mm=p.mm,p.follow=p.follow,n.add=n.add)
		}
	}
	return
}