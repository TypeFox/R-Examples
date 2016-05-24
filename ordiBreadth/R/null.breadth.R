null.breadth <-
function(dat,dist.method="jaccard",rep=100,quantiles=c(0.025,0.975),scaled=FALSE){
	hstrich<-sort(unique(apply(dat,1,sum)))
	nulls<-array(dim=c(length(hstrich),rep))
	for(i in 1:length(hstrich)){
		tosamp<-rep(0,dim(dat)[2]);tosamp[1:hstrich[i]]<-1
		cat("\n",i,"of",length(hstrich))
		for(j in 1:rep){
			nulls[i,j]<-hyp.ordi.breadth(dat,sample(tosamp,size=length(tosamp)))
					}
				}
			cat("\n")	
			u.q<-NA
			l.q<-NA	
			for(i in 1:dim(nulls)[1]){
				u.q[i] <- quantile(nulls[i,],quantiles[2])
				l.q[i] <- quantile(nulls[i,],quantiles[1])
			}
				if(scaled==TRUE){
					ug<-hyp.ordi.breadth(dat,rep(1,length(tosamp)),dist.method=dist.method) 
					u.q<-u.q/ug;l.q<-l.q/ug
							}
			richness<-hstrich				
			res<-cbind(richness,l.q,u.q)	
			return(res)
		}
