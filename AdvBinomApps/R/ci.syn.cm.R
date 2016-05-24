ci.syn.cm <- function(k,n,K,theta,alpha=0.1,p.target=1,tol=1e-10){
	
	if((alpha<=0)||(alpha>=1)||(p.target<=0)||(p.target>1)){
		return<-"Input error: check alpha and/or p.target!"
	}else{
		phi.cm<-phi.syn.cm(k,n,K,theta)
		if(is.character(phi.cm)){
			return<-phi.cm
		}else{
			if(phi.cm$prob[1]==1){
				p.hat<-qbeta(1-alpha,phi.cm$k+1,min(n)-phi.cm$k)
			}else{
				f.p<-function(p){as.numeric(phi.cm$prob%*%pbinom(phi.cm$k,min(n),p))-alpha}
     				p.hat<-uniroot(f.p,lower=0,upper=1,tol=tol)$root
      		}

			if(p.hat>p.target){	
				f.n<-function(n.add){
					n.new<-floor(n+n.add)
					phi.cm<-phi.syn.cm(k,n.new,K,theta)
					if(phi.cm$prob[1]==1){
						p<-qbeta(1-alpha,phi.cm$k+1,min(n.new)-phi.cm$k)
					}else{
						f.p<-function(p){phi.cm$prob%*%pbinom(phi.cm$k,min(n.new),p)-alpha}
     						p<-uniroot(f.p,lower=0,upper=1,tol=tol)$root
					}
					p-p.target
				}
				n.add<-ceiling(uniroot(f.n,lower=0,upper=100*max(n),extendInt="yes")$root)
			}else{
				n.add<-0
			}	
			return<-list(p.hat.cm=p.hat,n.add.cm=n.add)
		}	
	}
	return
}