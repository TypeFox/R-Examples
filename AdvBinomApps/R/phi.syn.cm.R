phi.syn.cm <- function(k,n,K,theta){

	r<-length(k)
	
	if(length(n)!=length(k)){return<-"Input error: dimensions of k and n do not match!"
	}else if(min(n-k)<0){return<-"Input error: k has to be <= n!"
	}else if(min(n)<=0){return<-"Input error: n has to be > 0!"
	}else if(!is.matrix(K)){return<-"Input error: K has to be a matrix!"
	}else if((r!=dim(K)[1])||(length(theta)!=dim(K)[2])){return<-"Input error: dimensions of k and K and/or theta and K do not match!"
	}else if((min(theta)<0)||(max(theta)>1)){return<-"Input error: theta has to be in [0,1]!"
	}else{

		K<-cbind(K,k-apply(K,1,sum))
		theta<-c(theta,0)

		if(min(K[,length(theta)])<0){
			return<-"Input error: arguments in k and K do not match!"
		}else{
			if(r>1){
				L<-as.matrix(c(0:k[r]))
				for(i in c((r-1):1)){
					L.new<-L
					if(k[i]>0){
						for(j in c(1:k[i])){
							L.new<-rbind(L.new,L)
						}
					}	
					L<-cbind(sort(rep(c(0:k[i]),dim(L.new)[1]/(k[i]+1))),L.new)
				}
	
				phi.L<-NULL
				for(i in c(1:dim(L)[1])){
					phi.syn.i<-phi.syn(L[i,],n)
		
					prob.i<-1
					for(j in c(1:r)){
						if(k[j]>0){
							prob.i<-prob.i*dgbinom(L[i,j],K[j,],1-theta)	
						}
					}

					phi.syn.i$prob<-phi.syn.i$prob*prob.i
					phi.L<-rbind(phi.L,phi.syn.i)
				}
	
				phi.cm<-NULL
				omega.r<-c(0:max(phi.L$k))
				for(j in omega.r){
					phi.cm<-c(phi.cm,sum(phi.L$prob[phi.L$k==j]))
				}

				omega.r<-omega.r[phi.cm>0]
				phi.cm<-phi.cm[phi.cm>0]

				return<-data.frame(k=omega.r,prob=phi.cm)

			}else{
				if(k==0){
					return<-data.frame(k=0,prob=1)
				}else{
					prob<-dgbinom(c(0:k),K[1,],1-theta)
					return<-data.frame(k=c(0:k)[prob>0],prob=prob[prob>0])
				}
			}
		}	
	}
	return
}
