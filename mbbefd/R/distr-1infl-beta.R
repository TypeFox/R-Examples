#d, p, q, r function for one-inflated beta distribution


doibeta <- function(x, shape1, shape2, p1, ncp=0, log=FALSE)
{
  doifun(x=x, dfun=dbeta, p1=p1, log=log, shape1=shape1, shape2=shape2, ncp=ncp)
}

poibeta <- function(q, shape1, shape2, p1, ncp=0, lower.tail = TRUE, log.p = FALSE)
{
  poifun(q=q, pfun=pbeta, p1=p1, lower.tail = lower.tail, log.p = log.p, shape1=shape1, shape2=shape2, ncp=ncp)
}

qoibeta <- function(p, shape1, shape2, p1, ncp=0, lower.tail = TRUE, log.p = FALSE)
{
  qoifun(p=p, qfun=qbeta, p1=p1, lower.tail = lower.tail, log.p = log.p, shape1=shape1, shape2=shape2, ncp=ncp)
}  

roibeta <- function(n, shape1, shape2, p1, ncp=0)
{
  roifun(n=n, rfun=rbeta, p1=p1, shape1=shape1, shape2=shape2, ncp=ncp)
}

ecoibeta <- function(x, shape1, shape2, p1, ncp=0)
{
  if(ncp != 0)
    stop("not yet implemented.")
  ecoifun(x=x, ecfun=ecbeta, mfun=mbeta, p1=p1, shape1=shape1, shape2=shape2)
}

moibeta <- function(order, shape1, shape2, p1, ncp=0)
{
  if(ncp != 0)
    stop("not yet implemented.")
  moifun(order=order, mfun=mbeta, p1=p1, shape1=shape1, shape2=shape2)
}


tloibeta <- function(shape1, shape2, p1, ncp=0)
{
  tloifun(p1=p1, shape1=shape1, shape2=shape2, ncp=ncp)
}
