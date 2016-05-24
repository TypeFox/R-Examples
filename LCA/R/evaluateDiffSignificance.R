evaluateDiffSignificance <- function(d,diff,PTLmodel){
        model_params <- predictPTLparams(d,PTLmodel)
	test.cdf <- ecdf(rep(seq(from=-1,to=1,by=0.01),round(100*unlist(lapply(seq(from=-1,to=1,by=0.01),dPTL,alpha=model_params$alpha,beta=model_params$beta,gamma=model_params$gamma)))))
        if(diff<0){
                out <- test.cdf(diff)
        }
        else{
                out <- 1-test.cdf(diff)
        }
        out
}
