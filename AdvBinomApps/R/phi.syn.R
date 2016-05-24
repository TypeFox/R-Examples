phi.syn <- function(k,n){
	
	if(length(n)!=length(k)){return<-"Input error: dimensions of k and n do not match!"
	}else if(min(n-k)<0){return<-"Input error: k has to be <= n!"
	}else if(min(n)<=0){return<-"Input error: n has to be > 0!"
	}else{
		r<-length(k)
		if(r>1){
			ord<-order(n,decreasing=TRUE)
			n<-n[ord]
			k<-k[ord]
			n.min<-min(n)
			omega.r<-c(max(k-(n-n.min)):min(sum(k),n.min))

			prob.k.kstar<-function(k,k.star){
				if((sum(k.star)<k)||(max(k.star)>k)){
					return<-0
				}else{
					if(r==2){
						return<-choose(n.min,k.star[1])*choose(n.min-k.star[1],k-k.star[1])*choose(k.star[1],k.star[1]+k.star[2]-k)/(choose(n.min,k.star[1])*choose(n.min,k.star[2]))
					}else{
						ord<-order(k.star,decreasing=TRUE)
						k.star<-k.star[ord]
						J.star<-as.matrix(c(0:k.star[r-1]))
						if(r>3){
							for(i in c((r-2):2)){
								J.star.new<-J.star
								if(k.star[i]!=0){	
									for(j in c(1:k.star[i])){
										J.star.new<-rbind(J.star.new,J.star)
									}
								}
								J.star<-cbind(sort(rep(c(0:k.star[i]),dim(J.star.new)[1]/(k.star[i]+1))),J.star.new)
							}
						}
						J.star<-cbind(J.star,k-apply(J.star,1,sum)-k.star[1])
						J.star<-J.star[(J.star[,r-1]>=0)&(J.star[,r-1]<=k.star[r-1]),]
						if(is.null(dim(J.star))){
							J.star<-matrix(J.star,1,r-1)	
						}
						f2apply<-function(rowK){
							choose(n.min,k.star[1])*prod(choose(n.min-cumsum(c(k.star[1],rowK[-(r-1)])),rowK)*choose(cumsum(c(k.star[1],rowK[-(r-1)])),k.star[-1]-rowK))
						}
						return<-sum(apply(J.star,1,f2apply))/prod(choose(n.min,k.star))
					}
				}
				return
			}

			K.star<-as.matrix(c(max(0,k[r-1]-(n[r-1]-n.min)):min(k[r-1],n.min)))
			if(r>2){
				for(i in c((r-2):1)){
					T.i<-c(max(0,k[i]-(n[i]-n.min)):min(k[i],n.min))
					K.star.new<-K.star
					if(length(T.i)>1){
						for(j in c(1:(length(T.i)-1))){
							K.star.new<-rbind(K.star.new,K.star)
						}
					}
					K.star<-cbind(sort(rep(T.i,dim(K.star.new)[1]/length(T.i))),K.star.new)
				}
			}
			K.star.mod<-cbind(K.star,rep(k[r],dim(K.star)[1]))
			phi<-NULL
			for(j in omega.r){
				phi.j<-0
				for(i in c(1:dim(K.star.mod)[1])){
					phi.j<-phi.j+prob.k.kstar(j,K.star.mod[i,])*prod(dhyper(K.star.mod[i,c(1:(r-1))],k[-r],n[-r]-k[-r],n.min))
				}
				phi<-c(phi,phi.j)
			}
			return<-data.frame(k=omega.r,prob=phi)
		}else{
			return<-data.frame(k=k,prob=1)	
		}
	}
	return
}
