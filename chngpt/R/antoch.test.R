antoch.test=function(formula, data, chngpt.var, plot.=FALSE) {
    
    DNAME = deparse(substitute(data))
    
    dat.1=data[order(data[[chngpt.var]]),]
    
    tmp=model.frame(formula, dat.1)
    y=tmp[[1]]
    
    fit=keepWarnings(glm(formula,  data=dat.1, family="binomial"))    # if glm gives a warning, use sandwich estimator to get variance est
    if(length(fit$warning)!=0) {
        return (NA)
    } else {
        fit=fit$value
    }
    
    n = nrow(dat.1)
    mu.hat = expit(predict(fit))
    
    T. = sapply (1:(n-1), function(k) {
        S.k.0 = sum((y-mu.hat)[1:k])
        V.hat = mean(mu.hat*(1-mu.hat)) * k * (n-k) / n # this V.hat is a simplified version of the actual variance if the formula has more than an intercept as predictors
        abs(S.k.0)/sqrt(V.hat)
    })    
    if(plot.) plot(T., type="b")
    
    # compute p-value
    T.max=max(T.)
    names(T.max)="Maximum statisticsZ"
    loglogn=log(log(n))
    T1=sqrt(2*loglogn)*T.max - 2*loglogn - 1/2*log(loglogn) + 1/2*log(pi)
    #p.value=exp(-2*exp(-T1)) 
    p.value=1-exp(-2*exp(-T1)) 
    
    res=list()
    res$statistic=T.max
    res$p.value=p.value
    res$k=which.max(T.)
    res$parameter=NULL
    res$conf.int=NULL
    res$estimate=NULL
    res$null.value=NULL
    res$alternative="two-sided"
    res$method="Antoch Change Point Test"
    res$data.name=DNAME
        
    class(res)=c("htest",class(res))
    res
}
