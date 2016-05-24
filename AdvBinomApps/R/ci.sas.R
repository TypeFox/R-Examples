ci.sas <- function(k, n, A.ref, A.follow, alpha=0.1, p.target=1, atol=1e-08){
	
	if((alpha<=0)||(alpha>=1)||(p.target<=0)||(p.target>1)){return<-"Input error: check alpha and/or p.target!"
	}else if((length(k)!=length(A.ref))||(length(A.ref)!=length(A.follow))){return<-"Input error: dimensions of k, A.ref and A.follow do not match!"
	}else if(n-sum(k)<0){return<-"Input error: total number of failures has to be <= n!"
	}else if(n<=0){return<-"Input error: n has to be > 0!"
	}else if(min(A.ref)<=0){return<-"Input error: A.ref has to be > 0!"
	}else if(min(A.follow)<0){return<-"Input error: A.follow has to be >= 0!"
	}else{

		r <- length(k)
		k.tot<-sum(k)
		A.ref.tot <- sum(A.ref)
		A.follow.tot <- sum(A.follow)	

		p<-qbeta(1-alpha,k.tot+1,n-k.tot)

		#CAS
		p.cas<-1-(1-p)^(A.ref/A.ref.tot)
		p.follow.cas<-1-(1-p)^(A.follow.tot/A.ref.tot)

		if(p.follow.cas>p.target){
			f.n<-function(n.new){pbinom(k.tot,floor(n.new),1-(1-p.target)^(A.ref.tot/A.follow.tot))-alpha}
			n.add.cas<-ceiling(uniroot(f.n,lower=n,upper=100*n,extendInt="yes")$root)-n
		}else{
			n.add.cas<-0
		}

		if(r>1){

			#delta
			delta.ij<-NULL
			for(i in c(1:(r-1))){
				for(j in c((i+1):r)){
					if(k[i]>=k[j]){
						if((k[j]==0)||(A.ref[i]<=A.ref[j])){
							delta.ij<-c(delta.ij,pbinom(k[i],n,p.cas[i])/pbinom(k[j],n,p.cas[j]))
						}else{
							delta.ij<-c(delta.ij,pbinom(k[j],n,p.cas[j])/pbinom(k[i],n,p.cas[i]))
						}
					}else{
						if((k[i]==0)||(A.ref[j]<=A.ref[i])){
							delta.ij<-c(delta.ij,pbinom(k[j],n,p.cas[j])/pbinom(k[i],n,p.cas[i]))
						}else{
							delta.ij<-c(delta.ij,pbinom(k[i],n,p.cas[i])/pbinom(k[j],n,p.cas[j]))
						}
					}
				}
			}
			delta<-max(delta.ij)

			#SAS
			f.sas<-function(p.sas){c((1-p)-prod(1-p.sas),pbinom(k[-r],n,p.sas[-r])-pbinom(k[-1],n,p.sas[-1]))}
			start.val<-1-(1-p)^((k+1)/sum(k+1))
			p.sas<-multiroot(f.sas,start.val,positive=TRUE,atol=atol,ctol=0,rtol=0)
			if(max(abs(p.sas$f.root))>atol){
				p.sas<-rep(NA,r)	
				p.follow.sas<-NA	
				n.add.sas<-NA
			}else{
				p.sas<-p.sas$root
				p.follow.sas<-1-prod((1-p.sas)^(A.follow/A.ref))

				if(p.follow.sas>p.target){
					f.n<-function(n.new){
						p<-qbeta(1-alpha,k.tot+1,floor(n.new)-k.tot)
						start.val<-1-(1-p)^((k+1)/sum(k+1))
						f.sas<-function(p.sas){c((1-p)-prod(1-p.sas),pbinom(k[-r],floor(n.new),p.sas[-r])-pbinom(k[-1],floor(n.new),p.sas[-1]))}
						p.sas<-multiroot(f.sas,start.val,positive=TRUE,atol=atol,ctol=0,rtol=0)$root
						p.follow.sas<-1-prod((1-p.sas)^(A.follow/A.ref))
						p.follow.sas-p.target
					}
					n.add.sas<-ceiling(uniroot(f.n,lower=n,upper=100*n,extendInt="yes")$root)-n
				}else{
					n.add.sas<-0
				}
			}

		}else{
			p.sas<-p.cas
			p.follow.sas<-p.follow.cas
			delta<-1
			n.add.sas<-n.add.cas		
		}
		return<-list(p.cas=p.cas, p.sas=p.sas, p.follow.cas=p.follow.cas, p.follow.sas=p.follow.sas, delta=delta, n.add.cas=n.add.cas,n.add.sas=n.add.sas)
	}
	return
}
