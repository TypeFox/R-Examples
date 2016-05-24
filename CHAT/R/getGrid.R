getGrid <- function(x0,y0,p,type=1,Nc=10,BAF='all',para=para){
	# (x0,y0) new coordination for origin
	# p purity for BAF
	# q1 contraction for logR on deletion side
	# q2 contraction for logR on amplification side
	# Nc maxmum copy number desired
	# keep track of BAF=0,1,2
    
	nn<-type
	if(para$model==1){
		n0<-log2(1+(nn-1)*p)+1
		nn<-1
	}
	else{
		n0<-log2(nn)+1
	}
	if(para$is.tri){
		n0<-log2(2+p)
		b0<-abs(1/(2+p)-0.5)
		nn<-1
	}
	gg<-c()
	if(BAF=='all'){
		for(i in 1:Nc){
			for(nc in 0:i){
				if(nc>i/2){next}
				ll<-(log2(p*i+(1-p)*2*nn)-n0)+y0
				bb<-abs((p*nc+(1-p)*nn)/(p*i+(1-p)*2*nn)-0.5)+x0
				if(para$is.tri)bb<-abs((p*nc+(1-p))/(p*i+(1-p)*2*nn)-0.5)+x0-b0
				gg<-rbind(gg,c(bb,ll))
			}
		}
	}
	else{
		for(i in 1:Nc){
			nc<-i-BAF
			if(nc<BAF){next}
			ll<-(log2(p*i+(1-p)*2*nn)-n0)+y0
			bb<-abs((p*nc+(1-p)*nn)/(p*i+(1-p)*2*nn)-0.5)+x0
			if(para$is.tri)bb<-abs((p*nc+(1-p))/(p*i+(1-p)*2*nn)-0.5)+x0-b0
			gg<-rbind(gg,c(bb,ll))
		}
	}
	return(gg)
}