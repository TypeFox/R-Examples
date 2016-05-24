# one-inflated shifted truncated Pareto distribution

doistpareto <- function(x, a, p1, log=FALSE)
{
  doifun(x=x, dfun=dstpareto, p1=p1, log=log, a=a)
}

poistpareto <- function(q, a, p1, lower.tail = TRUE, log.p = FALSE)
{
  poifun(q=q, pfun=pstpareto, p1=p1, lower.tail = lower.tail, log.p = log.p, a=a)
}

qoistpareto <- function(p, a, p1, lower.tail = TRUE, log.p = FALSE)
{
  qoifun(p=p, qfun=qstpareto, p1=p1, lower.tail = lower.tail, log.p = log.p, a=a)
}  

roistpareto <- function(n, a, p1)
{
  roifun(n=n, rfun=rstpareto, p1=p1, a=a)
}

ecoistpareto <- function(x, a, p1)
{
  ecoifun(x=x, ecfun=ecstpareto, mfun=mstpareto, p1=p1, a=a)
}

moistpareto <- function(order, a, p1)
{
  moifun(order=order, mfun=mstpareto, p1=p1, a=a)
}

tloistpareto <- function(a, p1)
{
  tloifun(p1=p1, a=a)
}

