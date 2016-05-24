#d, p, q, r function for one-inflated uniform distribution


doiunif <- function(x, p1, log=FALSE)
{
  doifun(x=x, dfun=dunif, p1=p1, log=log, min=0, max=1)
}

poiunif <- function(q, p1, lower.tail = TRUE, log.p = FALSE)
{
  poifun(q=q, pfun=punif, p1=p1, lower.tail = lower.tail, log.p = log.p, min=0, max=1)
}

qoiunif <- function(p, p1, lower.tail = TRUE, log.p = FALSE)
{
  qoifun(p=p, qfun=qunif, p1=p1, lower.tail = lower.tail, log.p = log.p, min=0, max=1)
}  

roiunif <- function(n, p1)
{
  roifun(n=n, rfun=runif, p1=p1, min=0, max=1)
}


ecoiunif <- function(x, p1)
{
  ecoifun(x=x, ecfun=ecunif, mfun=munif, p1=p1, min=0, max=1)
}

moiunif <- function(order, p1)
{
  moifun(order=order, mfun=munif, p1=p1, min=0, max=1)
}

tloiunif <- function(p1)
{
  tloifun(p1=p1)
}
