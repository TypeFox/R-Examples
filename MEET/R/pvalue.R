pvalue <-function (a,dist){
	require("KernSmooth")
	require("Hmisc")
    
    
if(sd(dist)==0){
	if(a>mean(dist)){
		pvalue <- 0
	}else{
		pvalue <- 1	
	}
}else{
	h <- dpik(as.vector(dist))
	est <- bkde(dist, bandwidth=h,range.x=extendrange(r=range(a,dist), f = 0.05), gridsize=1001)
	pvalue <- trap.rule(est$x[est$x>a], est$y[est$x>a])/trap.rule(est$x, est$y)}
pvalue <- length(dist[dist>a])/length(dist)
pvalue
}

