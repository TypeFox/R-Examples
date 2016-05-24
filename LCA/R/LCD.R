LCD <- function(x1,x2,seed.row,PTLmodel){
	d <- sqrt(sum((x1-x2)^2,na.rm=TRUE))

        model_params <- predictPTLparams(d,PTLmodel)
        test.cdf <- ecdf(rep(seq(from=-1,to=1,by=0.01),round(100*unlist(lapply(seq(from=-1,to=1,by=0.01),dPTL,alpha=model_params$alpha,beta=model_params$beta,gamma=model_params$gamma)))))

	mapply(function(x,y){max(c(2.2e-16,test.cdf(max(x,y)[1])-test.cdf(min(x,y)[1])))},x=x2-x1,y=rep(x2[seed.row]-x1[seed.row],length(x1)))
}
