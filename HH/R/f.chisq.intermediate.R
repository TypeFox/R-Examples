## intermediate functions to simplify writing for both R and S-Plus

dchisq.intermediate <- function(x, df, ncp=0, log=FALSE)
  if.R(r=dchisq(x=x, df=df, ncp=ncp, log=log),
       s=dchisq(x=x, df=df, log=log))

pchisq.intermediate <- function(q, df, ncp=0, lower.tail=TRUE, log.p=FALSE)
  if.R(r=pchisq(q=q, df=df, ncp=ncp, lower.tail=lower.tail, log.p=log.p),
       s=pchisq(q=q, df=df, ncp=ncp))

qchisq.intermediate <- function(p, df, ncp=0, lower.tail=TRUE, log.p=FALSE)
  if.R(r=qchisq(p=p, df=df, ncp=ncp, lower.tail=lower.tail, log.p=log.p),
       s=qchisq(p=p, df=df))


df.intermediate <- function(x, df1, df2, ncp=0, log=FALSE)
  if.R(r=df(x=x, df1=df1, df2=df2, ncp=ncp, log=log),
       s=df(x=x, df1=df1, df2=df2, log=log))

pf.intermediate <- function(q, df1, df2, ncp=0, lower.tail=TRUE, log.p=FALSE)
  if.R(r=pf(q=q, df1=df1, df2=df2, ncp=ncp, lower.tail=lower.tail, log.p=log.p),
       s=pf(q=q, df1=df1, df2=df2, ncp=ncp))

qf.intermediate <- function(p, df1, df2, ncp=0, lower.tail=TRUE, log.p=FALSE)
  if.R(r=qf(p=p, df1=df1, df2=df2, ncp=ncp, lower.tail=lower.tail, log.p=log.p),
       s=qf(p=p, df1=df1, df2=df2))

## source("~/HH-R.package/HH/R/f.chisq.intermediate.R")

## ## R
## df(x, df1, df2, ncp, log = FALSE)
## pf(q, df1, df2, ncp, lower.tail = TRUE, log.p = FALSE)
## qf(p, df1, df2, ncp, lower.tail = TRUE, log.p = FALSE)

## dchisq(x, df, ncp=0, log = FALSE)
## pchisq(q, df, ncp=0, lower.tail = TRUE, log.p = FALSE)
## qchisq(p, df, ncp=0, lower.tail = TRUE, log.p = FALSE)

## ## S-Plus
## df(x, df1, df2, log=F) 
## pf(q, df1, df2, ncp=0) 
## qf(p, df1, df2) 

## dchisq(x, df, log=F) 
## pchisq(q, df, ncp=0) 
## qchisq(p, df) 
