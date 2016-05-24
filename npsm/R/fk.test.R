fk.test <-
function(x,y,alternative = c("two.sided", "less", "greater"),conf.level = 0.95){
    fkscores = new("scores",phi=function(u)
    {
        (qnorm((u+1)/2))^2 - 1
    }
    ,Dphi=function(u)
    {
        qnorm((u+1)/2)/dnorm(qnorm((u+1)/2))
    })
#    myscores = fkscores
#
#      Test
#
    zed = c(abs(x-median(x)),abs(y-median(y)))
    n1 = length(x)
    n2 = length(y)
    n = n1 + n2
    cvec  = c(rep(0,n1),rep(1,n2))
    rz = rank(zed)/(n+1)
    v = getScores(fkscores,rz)
    num = sum(cvec*v)
    sigphi = sqrt(((n1*n2)/n)*var(v))
    ts = num/sigphi
    see = alternative[1]
    if(see == "two.sided"){pval = 2*(1 - pnorm(abs(ts)))}
    if(see == "greater"){pval = 1 - pnorm(ts)}
    if(see == "less"){pval = pnorm(ts)}

#
#      Estimation and CI
#
    xs = abs(x - median(x))
    xs = xs[xs!=0]
    xstarl = log(xs)
    ys = abs(y - median(y))
    ys = ys[ys!=0]
    ystarl = log(ys)
    n1s = length(xs)
    n2s = length(ys)
    ns = n1s + n2s
    tc = abs(qt(((1-conf.level)/2),ns-2))
    cvec  = c(rep(0,n1s),rep(1,n2s))
    zed = c(xstarl,ystarl)
    fitz = rfit(zed~cvec,scores=fkscores)
    sumf = summary(fitz)
    delta = coef(sumf)[2,1]
    se = coef(sumf)[2,2]
    lb = delta - tc*se
    ub = delta + tc*se
    eta = exp(delta) 
    ci = c(exp(lb),exp(ub))
    res<-list(estimate=eta,conf.int=ci,statistic=ts,p.value=pval,conf.level=conf.level)
    res$call<-match.call()
    class(res)<-list('rank.test')
    res
}
