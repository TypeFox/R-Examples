phi.mult.ref <- function(k,n,A.ref,prec=2,tailcut=1e-08){	

	if((length(n)!=length(k))||(length(k)!=length(A.ref))){return<-"Input error: dimensions of k, n and A.ref do not match!"
	}else if(min(n-k)<0){return<-"Input error: k has to be <= n!"
	}else if(min(n)<=0){return<-"Input error: n has to be > 0!"
	}else if((as.character(prec)!=as.character(round(prec)))||(prec<0)){return<-"Input error: prec has to be a positive integer!"	
	}else if(tailcut<0){return<-"Input error: tailcut has to be >= 0!"
	}else if(min(A.ref)<10^(-prec)){return<-"Input error: A.ref has to be >= 10^-prec!"
	}else{

		A.prec<-round(A.ref,prec)	
		A.prec[A.prec-A.ref>0]<-A.prec[A.prec-A.ref>0]-10^(-prec)
		A.ref<-A.prec

		prob.uj.kj.gcd<-function(uj,kj.gcd,nj.gcd,mj.gcd){
			if(uj>kj.gcd){
				return<-0
			}else if((uj==0)&(kj.gcd>0)){
				return<-0
			}else{
				if(uj>1){
					Y<-as.matrix(c(1:(min(kj.gcd-uj+1,mj.gcd))))
					if(uj>2){
						for(j in c(2:(uj-1))){
							Y.new<-Y
							if(kj.gcd>uj){
								for(i in c(2:min(kj.gcd-uj+1,mj.gcd))){
									Y.new<-rbind(Y.new,Y)
								}
							}
							Y<-cbind(sort(rep(c(1:min(kj.gcd-uj+1,mj.gcd)),dim(Y.new)[1]/min(kj.gcd-uj+1,mj.gcd))),Y.new)
						}
					}
					Y<-cbind(Y,kj.gcd-apply(Y,1,sum))
					Y<-Y[(Y[,uj]>0)&(Y[,uj]<=mj.gcd),]
				}else{
					Y<-kj.gcd
				}

				if(is.null(dim(Y))){
					Y<-matrix(Y,1,length(Y))
				}

				f.prob<-function(rowY){
					prod(dhyper(rowY,kj.gcd-c(0,cumsum(rowY)[-uj]),nj.gcd-c(0:(uj-1))*mj.gcd-(kj.gcd-c(0,cumsum(rowY)[-uj])),mj.gcd))
				}

				return<-choose(nj.gcd/mj.gcd,uj)*sum(apply(Y,1,f.prob))
			}
			return
		}

		phi.j<-function(kj,nj,Aj,A.gcd){
			mj.gcd<-Aj/A.gcd
			nj.gcd<-nj*mj.gcd
			Tj.gcd<-c(kj:(kj*mj.gcd))
			cum.probs<-1
			if(length(Tj.gcd)>1){	
				for(kj.gcd in Tj.gcd[-1]){
					supp.uj<-c(ceiling(kj.gcd/mj.gcd):kj)
					uj<-matrix(supp.uj,length(supp.uj),1)
					cum.prob.kj.gcd<-sum(apply(uj,1,prob.uj.kj.gcd,kj.gcd,nj.gcd,mj.gcd))
					if(cum.prob.kj.gcd<tailcut){break}else{cum.probs<-c(cum.probs,cum.prob.kj.gcd)}
				}
				if(length(cum.probs)>1){
					phi<-cum.probs-c(cum.probs[2:length(cum.probs)],0)
				}else{
					phi<-NULL
				}
			}else{
				phi<-cum.probs
			}
			if(is.null(phi)){
				return<-"Input error: tailcut too large!"
			}else{
				return<-data.frame(kj.gcd=Tj.gcd[1:length(phi)],prob=phi)
			}
			return
		}

		r<-length(k)
		A.gcd<-gcd.mult.ref(A.ref,prec)
		phi1<-phi.j(k[1],n[1],A.ref[1],A.gcd)
	
		if(r>1){
			if(is.character(phi1)){
				return<-phi1
			}else{
				for(j in c(2:r)){
					phi2<-phi.j(k[j],n[j],A.ref[j],A.gcd)
					if(is.character(phi2)){
						phi1<-phi2
						break
					}else{
						T.gcd<-c((phi1$kj.gcd[1]+phi2$kj.gcd[1]):(max(phi1$kj.gcd)+max(phi2$kj.gcd)))
						prob.k.gcd<-NULL
						for(k.gcd in T.gcd){
							k1.gcd<-c(max(phi1$kj.gcd[1],k.gcd-max(phi2$kj.gcd)):min(k.gcd-phi2$kj.gcd[1],max(phi1$kj.gcd)))
							k2.gcd<-k.gcd-k1.gcd
							prob.k.gcd<-c(prob.k.gcd,sum(phi1$prob[k1.gcd-(phi1$kj.gcd[1]-1)]*phi2$prob[k2.gcd-(phi2$kj.gcd[1]-1)]))
						}
						phi1<-data.frame(kj.gcd=T.gcd,prob=prob.k.gcd)	
					}
				}	
			}
		}	
		if(is.character(phi1)){
			return<-phi1
		}else{
			dimnames(phi1)[[2]][1]<-"k.gcd"
			return<-list(phi=phi1,A.gcd=A.gcd)
		}
	}
	return
}
