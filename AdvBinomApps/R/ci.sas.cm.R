ci.sas.cm <- function(k, n, A.ref, A.follow, K, theta, alpha=0.1, p.target=1, atol=1e-08){
	
	r<-length(k)

	if((alpha<=0)||(alpha>=1)||(p.target<=0)||(p.target>1)){return<-"Input error: check alpha and/or p.target!"
	}else if((length(k)!=length(A.ref))||(length(A.ref)!=length(A.follow))){return<-"Input error: dimensions of k, n, A.ref and A.follow do not match!"
	}else if(n-sum(k)<0){return<-"Input error: total number of failures has to be <= n!"
	}else if(n<=0){return<-"Input error: n has to be > 0!"
	}else if(min(A.ref)<=0){return<-"Input error: A.ref has to be > 0!"
	}else if(min(A.follow)<0){return<-"Input error: A.follow has to be >= 0!"
	}else if(!is.matrix(K)){return<-"Input error: K has to be a matrix!"
	}else if((r!=dim(K)[1])||(length(theta)!=dim(K)[2])){return<-"Input error: dimensions of k and K and/or theta and K do not match!"
	}else if((min(theta)<0)||(max(theta)>1)){return<-"Input error: theta has to be in [0,1]!"
	}else{

		k.tot<-sum(k)

		if(k.tot==0){
			return<-ci.sas(k,n,A.ref,A.follow,alpha,p.target,atol)
		}else{
	
			A.ref.tot<-sum(A.ref)
			A.follow.tot<-sum(A.follow)	

			k.with.cm<-apply(K,1,sum)
			k.without.cm<-k-k.with.cm
			K<-cbind(K,k.without.cm)
			theta<-c(theta,0)

			if(min(K[,length(theta)])<0){
				return<-"Input error: arguments in k and K do not match!"
			}else{
				p.cm<-cm.clopper.pearson.ci(n,apply(K,2,sum),theta,alpha,uniroot.tol=1e-12)$Upper.limit

				#CAS
				p.cas.cm<-1-(1-p.cm)^(A.ref/A.ref.tot)
				p.follow.cas.cm<-1-(1-p.cm)^(A.follow.tot/A.ref.tot)

				if(p.follow.cas.cm>p.target){
					f.n<-function(n.new){
						p.cm<-cm.clopper.pearson.ci(floor(n.new),apply(K,2,sum),theta,alpha,uniroot.tol=1e-12)$Upper.limit
						(1-(1-p.cm)^(A.follow.tot/A.ref.tot))-p.target
					}
					n.add.cas.cm<-ceiling(uniroot(f.n,lower=n,upper=100*n,extendInt="yes")$root)-n
				}else{
					n.add.cas.cm<-0
				}

				if(r>1){

					#delta
					delta.ij.cm<-NULL
					for(i in c(1:(r-1))){
						for(j in c((i+1):r)){

							if(k[i]>0){xi.i<-dgbinom(c(0:k[i]),K[i,],1-theta)}else{xi.i<-1}
							if(k[j]>0){xi.j<-dgbinom(c(0:k[j]),K[j,],1-theta)}else{xi.j<-1}
			
							if(k[i]>=k[j]){
								if((k[j]==0)||(A.ref[i]<=A.ref[j])){
									delta.ij.cm<-c(delta.ij.cm,xi.i%*%pbinom(c(0:k[i]),n,p.cas.cm[i])/xi.j%*%pbinom(c(0:k[j]),n,p.cas.cm[j]))
								}else{
									delta.ij.cm<-c(delta.ij.cm,xi.j%*%pbinom(c(0:k[j]),n,p.cas.cm[j])/xi.i%*%pbinom(c(0:k[i]),n,p.cas.cm[i]))
								}
							}else{
								if((k[i]==0)||(A.ref[j]<=A.ref[i])){
									delta.ij.cm<-c(delta.ij.cm,xi.j%*%pbinom(c(0:k[j]),n,p.cas.cm[j])/xi.i%*%pbinom(c(0:k[i]),n,p.cas.cm[i]))
								}else{
									delta.ij.cm<-c(delta.ij.cm,xi.i%*%pbinom(c(0:k[i]),n,p.cas.cm[i])/xi.j%*%pbinom(c(0:k[j]),n,p.cas.cm[j]))
								}
							}
						}
					}
					delta.cm<-max(delta.ij.cm)

					#SAS
					f.sas.cm<-function(p.sas.cm){
						eq<-numeric(r)
						eq[1]<-(1-p.cm)-prod(1-p.sas.cm)
						for(j in c(2:r)){
							if(k[j-1]>0){xi.jminus1<-dgbinom(c(0:k[j-1]),K[j-1,],1-theta)}else{xi.jminus1<-1}
							if(k[j]>0){xi.j<-dgbinom(c(0:k[j]),K[j,],1-theta)}else{xi.j<-1}
							eq[j]<-xi.jminus1%*%pbinom(c(0:k[j-1]),n,p.sas.cm[j-1])-xi.j%*%pbinom(c(0:k[j]),n,p.sas.cm[j])
						}
						eq
					}

					#Computation of appropriate starting values for multiroot-function
					support1<-c(0:k.with.cm[1])+k.without.cm[1]
					if(k.with.cm[1]>0){xi1<-dgbinom(c(0:k.with.cm[1]),K[1,-length(theta)],1-theta[-length(theta)])}else{xi1<-1}
					S<-as.matrix(support1[xi1>1e-05])
					Xi<-as.matrix(xi1[xi1>1e-05])
					for(i in c(2:r)){
						S.new<-NULL
						Xi.new<-NULL
						if(k.with.cm[i]>0){xi.i<-dgbinom(c(0:k.with.cm[i]),K[i,-length(theta)],1-theta[-length(theta)])}else{xi.i<-1}
						support.i<-(c(0:k.with.cm[i])+k.without.cm[i])[xi.i>1e-05]
						xi.i<-xi.i[xi.i>1e-05]
						for(j in support.i){
							S.new<-rbind(S.new,cbind(S,j))
							Xi.new<-rbind(Xi.new,cbind(Xi,xi.i[j-support.i[1]+1]))	
						}
						S<-S.new
						Xi<-Xi.new
					}
					xi.probs<-apply(Xi,1,prod)
					proportions<-(S+1)/(apply(S,1,sum)+r)
					M<-1-(1-p.cm)^proportions	
					start.val<-as.numeric(apply(M*xi.probs,2,sum))		
	
					p.sas.cm<-multiroot(f.sas.cm,start.val,positive=TRUE,atol=atol,ctol=0,rtol=0)
					if(max(abs(p.sas.cm$f.root))>atol){
						p.sas.cm<-rep(NA,r)	
						p.follow.sas.cm<-NA	
						n.add.sas.cm<-NA
					}else{
						p.sas.cm<-p.sas.cm$root
						p.follow.sas.cm<-1-prod((1-p.sas.cm)^(A.follow/A.ref))

						if(p.follow.sas.cm>p.target){
							f.n<-function(n.new){
								p.cm<-cm.clopper.pearson.ci(floor(n.new),apply(K,2,sum),theta,alpha,uniroot.tol=1e-12)$Upper.limit
								M<-1-(1-p.cm)^proportions	
								start.val<-as.numeric(apply(M*xi.probs,2,sum))	
								f.sas.cm<-function(p.sas.cm){
									eq<-numeric(r)
									eq[1]<-(1-p.cm)-prod(1-p.sas.cm)
									for(j in c(2:r)){
										if(k[j-1]>0){xi.jminus1<-dgbinom(c(0:k[j-1]),K[j-1,],1-theta)}else{xi.jminus1<-1}
										if(k[j]>0){xi.j<-dgbinom(c(0:k[j]),K[j,],1-theta)}else{xi.j<-1}
										eq[j]<-xi.jminus1%*%pbinom(c(0:k[j-1]),floor(n.new),p.sas.cm[j-1])-xi.j%*%pbinom(c(0:k[j]),floor(n.new),p.sas.cm[j])
									}
									eq
								}
								p.sas.cm<-multiroot(f.sas.cm,start.val,positive=TRUE,atol=atol,ctol=0,rtol=0)$root
								p.follow.sas.cm<-1-prod((1-p.sas.cm)^(A.follow/A.ref))
								p.follow.sas.cm-p.target
							}
							n.add.sas.cm<-ceiling(uniroot(f.n,lower=n,upper=100*n,extendInt="yes")$root)-n
						}else{
							n.add.sas.cm<-0
						}
					}
				}else{
					p.sas.cm<-p.cas.cm
					p.follow.sas.cm<-p.follow.cas.cm
					delta.cm<-1
					n.add.sas.cm<-n.add.cas.cm		
				}
				return<-list(p.cas.cm=p.cas.cm, p.sas.cm=p.sas.cm, p.follow.cas.cm=p.follow.cas.cm, p.follow.sas.cm=p.follow.sas.cm, delta.cm=delta.cm, n.add.cas.cm=n.add.cas.cm,n.add.sas.cm=n.add.sas.cm)
			}
		}
	}
	return
}
