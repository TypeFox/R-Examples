ci.syn <- function(k,n,alpha=0.1,p.target=1,tol=1e-10){

	if((alpha<=0)||(alpha>=1)||(p.target<=0)||(p.target>1)){
		return<-"Input error: check alpha and/or p.target!"
	}else{
		phi<-phi.syn(k,n)
		if(is.character(phi)){
			return<-phi
		}else{
			if(phi$prob[1]==1){
				p.hat<-qbeta(1-alpha,phi$k+1,min(n)-phi$k)
			}else{
      			f.p<-function(p){as.numeric(phi$prob%*%pbinom(phi$k,min(n),p))-alpha}
     				p.hat<-uniroot(f.p,lower=0,upper=1,tol=tol)$root
			}       
			if(p.hat>p.target){	
				f.n<-function(n.add){
					n.new<-floor(n+n.add)
					phi<-phi.syn(k,n.new)
					if(phi$prob[1]==1){
						p<-qbeta(1-alpha,phi$k+1,min(n.new)-phi$k)
					}else{
						f.p<-function(p){phi$prob%*%pbinom(phi$k,min(n.new),p)-alpha}
     						p<-uniroot(f.p,lower=0,upper=1,tol=tol)$root
					}
					p-p.target
				}
				n.add<-ceiling(uniroot(f.n,lower=0,upper=100*max(n),extendInt="yes")$root)
			}else{
				n.add<-0
			}
			return<-list(p.hat=p.hat,n.add=n.add)
		}
	}
	return
}
