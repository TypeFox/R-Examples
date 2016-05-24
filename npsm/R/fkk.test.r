fkk.test = function(y,ind,conf.level = 0.95){
    fkscores = new("scores",phi=function(u)
    {
        (qnorm((u+1)/2))^2 - 1
    }
    ,Dphi=function(u)
    {
        qnorm((u+1)/2)/dnorm(qnorm((u+1)/2))
    })
    myscores = fkscores
#
#      Test
#
    xmat = model.matrix(~as.factor(ind) - 1)
    uind <- unique(ind)
    k <- length(uind)
    yu = c(0)
    yua = c(0)
    iu = c(0)
    ni = c(0)
    for(j in 1:k){
       vec = y[xmat[,j]==1]
       yu = c(yu,vec)
       yua = c(yua,vec - median(vec))
       iu = c(iu,rep(j,length(vec)))
    }
    np1=length(yu)
    yu = yu[2:np1]
    yua = yua[2:np1]
    iu = iu[2:np1]
    regdata = cbind(yu,iu)
    n = np1 - 1
    rua = rank(abs(yua))/(n+1)
    v = getScores(myscores,(1:n)/(n+1))
    sv = sum(v^2)
    yuav = getScores(myscores,rua)
    ts = 0
    nc = c(0)
    for(j in 1:k){
       vec = yuav[xmat[,j]==1]
       ni = length(vec)
       nc = c(nc,ni)
       ts = ts + sum(vec)^2/ni
    }
    ts = ((n-1)/sv)*ts
    pval = 1 - pchisq(ts,k-1)

    nc  = nc[2:(k+1)]
    yua = c(0)
    iu = c(0)
    ni = c(0)
    for(j in 1:k){
       vec = y[xmat[,j]==1]
       vec = vec - median(vec)
       vec = vec[vec[]!=0]
       yua = c(yua,vec)
       iu = c(iu,rep(j,length(vec)))
    }
    np1=length(yua)
    yua = yua[2:np1]
    iu = iu[2:np1]
    n = np1 - 1
    zed = log(abs(yua))


    xm = model.matrix(~as.factor(iu))
    xmat = xm[,2:k]
    fitz = rfit(zed~xmat,scores=myscores)
    sumf = summary(fitz)
    delta = coef(sumf)[2:k,1]
    se = coef(sumf)[2:k,2]
    tc = abs(qt(((1-conf.level)/2),n-k))
    lb = delta - tc*se
    ub = delta + tc*se
    eta = exp(delta)
    ci = cbind(exp(lb),exp(ub))

    cwts = rep(1,nc[1])
    for(i in 2:k){
       cwts = c(cwts,rep(1/eta[i-1],nc[i]))
    }

    res<-list(statistic=ts,p.value=pval,estimate=eta,conf.int=ci,cwts=cwts,conf.level=conf.level)
    class(res)<-c('fkk.test')
    res
}
