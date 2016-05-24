#d, p, q, r function for one-inflated beta distribution


doigbeta <- function(x, shape0, shape1, shape2, p1, log=FALSE)
{
  doifun(x=x, dfun=dgbeta, p1=p1, log=log, shape1=shape1, shape2=shape2, shape0=shape0)
}

poigbeta <- function(q, shape0, shape1, shape2, p1, lower.tail = TRUE, log.p = FALSE)
{
  poifun(q=q, pfun=pgbeta, p1=p1, lower.tail=lower.tail, log.p=log.p, shape1=shape1, shape2=shape2, shape0=shape0)
}

qoigbeta <- function(p, shape0, shape1, shape2, p1, lower.tail = TRUE, log.p = FALSE)
{
  qoifun(p=p, qfun=qgbeta, p1=p1, lower.tail=lower.tail, log.p=log.p, shape1=shape1, shape2=shape2, shape0=shape0)
}  

roigbeta <- function(n, shape0, shape1, shape2, p1)
{
  roifun(n=n, rfun=rgbeta, p1=p1, shape1=shape1, shape2=shape2, shape0=shape0)
}

ecoigbeta <- function(x, shape0, shape1, shape2, p1)
{
  ecoifun(x=x, ecfun=ecgbeta, mfun=mgbeta, p1=p1, shape1=shape1, shape2=shape2, shape0=shape0)
}

moigbeta <- function(order, shape0, shape1, shape2, p1)
{
  moifun(order=order, mfun=mgbeta, p1=p1, shape1=shape1, shape2=shape2, shape0=shape0)
}

tloigbeta <- function(shape0, shape1, shape2, p1)
{
  tloifun(p1=p1, shape1=shape1, shape2=shape2, shape0=shape0)
}


