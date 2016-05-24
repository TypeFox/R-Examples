ci<-function(Beta,Beta_bootstrap,alpha,type=c("basic","quantile")){
# a=0 is corresponding to BC confidence interval	
	p<-dim(Beta_bootstrap)[2]
	times<-dim(Beta_bootstrap)[1]
	interval<-matrix(0,2,p)
	if(type=="basic"){
		bound.percentile<-apply(Beta_bootstrap,2,function(u){ quantile(u,prob=c(1-alpha/2,alpha/2)) })
		interval[1,]<-2*Beta-bound.percentile[1,]
		interval[2,]<-2*Beta-bound.percentile[2,]	
	}
	if(type=="quantile"){
		bound.percentile<-apply(Beta_bootstrap,2,function(u){ quantile(u,prob=c(alpha/2,1-alpha/2)) })
		interval[1,]<-bound.percentile[1,]
		interval[2,]<-bound.percentile[2,]		
	}

	return(interval)
}



