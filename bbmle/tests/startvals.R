library(bbmle)

## copied from emdbook
dbetabinom <- function (x, prob, size, theta, shape1, shape2, log = FALSE) 
{
    if (missing(prob) && !missing(shape1) && !missing(shape2)) {
        prob = shape1/(shape1 + shape2)
        theta = shape1 + shape2
    }
    v <- lchoose(size, x) - lbeta(theta * (1 - prob), theta * 
        prob) + lbeta(size - x + theta * (1 - prob), x + theta * 
        prob)
    if (log) 
        v
    else exp(v)
}

ss <- data.frame(taken=c(0,1,2,5),available=c(5,5,5,5),
                 dist=rep(1,4))

SP.bb=mle2(taken~dbetabinom(prob,theta,size=available),
  start=list(prob=0.5,theta=1),data=ss)
SP.bb.dist=mle2(taken~dbetabinom(prob,size=available,theta),
  parameters=list(prob~dist-1,theta~dist-1),
  start=as.list(coef(SP.bb)),data=ss)

SP.bb.dist2=mle2(taken~dbetabinom(prob,size=available,theta),
  parameters=list(prob~dist - 1,theta~dist - 1),
  start=as.list(coef(SP.bb)),data=ss)

