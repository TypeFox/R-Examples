################################################################################
################################################################################
Astar2<-function (d,n)
{
#  d <- as.matrix(d)
  m <- rowMeans(d)
  M <- mean(d)
  a <- sweep(d, 1, m)
  b <- sweep(a, 2, m)
  A <- b + M
  A <- A - d/n
  diag(A) <- m - M
  (n/(n - 1)) * A
}
################################################################################
# Description
# Computes distance correlation for functional data.
#
# Arguments
# D1: distances of 1st sample
# D2: distances of 2nd sample 
# Returns the sample distance correlation
dcor.dist=function(D1,D2){
       m1row=rowMeans(D1)
       m1col=colMeans(D1)
       m2row=rowMeans(D2)
       m2col=colMeans(D2)
       n=nrow(D1)
       ones=rep(1,n)
       pD1=D1-outer(ones,m1row)-outer(m1col,ones)+mean(D1)
       pD2=D2-outer(ones,m2row)-outer(m2col,ones)+mean(D2)
       res=sqrt(sum(pD1*pD2))/n
       res1=sqrt(sum(pD1*pD1))/n        
       res2=sqrt(sum(pD2*pD2))/n        
       out<-res/sqrt(res1*res2)
return(out)
}
################################################################################
################################################################################
# Description
# Computes the bias corrected distance correlation for functional data.
#
# Arguments
# D: distances of 1st sample
# D2: distances of 2nd sample 
# n:  sample dimension, by default nrow(D1)
# Returns the bias corrected dcor statistic
bcdcor.dist=function(D1,D2,n){
        m1row=rowMeans(D1)
        m1col=colMeans(D2)
        if (missing(n)) n<-nrow(D1)
        AA <- Astar2(D1,n)
        BB <- Astar2(D2,n)
        res <- sum(AA * BB) - (n/(n - 2)) * sum(diag(AA * BB))
        res1<- sum(AA * AA) - (n/(n - 2)) * sum(diag(AA * AA))
        res2<- sum(BB * BB) - (n/(n - 2)) * sum(diag(BB * BB))
        out<-res/sqrt(res1*res2)          
return(out)
}
################################################################################
################################################################################
dcor.test<-function (D1, D2,n)
{
  dname <- paste(deparse(substitute(D1)), "and", deparse(substitute(D2)))
  if (missing(n)) n <- nrow(D1)
  R<-bcdcor.dist(D1=D1,D2=D2,n=n) 
  M <- n * (n - 3)/2
  df <- M - 1
  names(df) <- "df"
  tstat <- sqrt(M - 1) * R/sqrt(1 - R^2)
  names(tstat) <- "T"
  names(R) <- "Bias corrected dcor"
  pval <- 1 - pt(tstat, df = df)
  method <- "dcor t-test of independence"
  rval <- list(statistic = tstat, parameter = df, p.value = pval,
               estimate = R, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}
################################################################################
################################################################################
dcor.xy<-function (x, y,test=TRUE,metric.x,metric.y,par.metric.x,par.metric.y,n)
{
  dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  if (missing(n)) n <- nrow(x)
  isfdata1<-is.fdata(x)
  true.dist<-FALSE
  if (missing(metric.x)){
    if (isfdata1) metric.x="metric.lp"
    else metric.x="metric.dist"
   }
  if (missing(par.metric.x)) par.metric.x<-list()
  if (missing(n)) n <- nrow(x)
  isfdata2<-is.fdata(y)
  if (isfdata1)         par.metric.x$fdata1<-x
  else                  par.metric.x$x<-as.matrix(x)
  D1=do.call(metric.x,par.metric.x)
  if (missing(metric.y)){
          if (isfdata2) metric.y="metric.lp"
          else metric.y="metric.dist"
  }
  if (missing(par.metric.y)) par.metric.y<-list()
  if (isfdata2)        par.metric.y$fdata1<-y
  else                 par.metric.y$x<-as.matrix(y)
  D2=do.call(metric.y,par.metric.y)
  D1<-as.matrix(D1)
  D2<-as.matrix(D2)
  if (test) {
   rval<-dcor.test(D1, D2,n)
   rval$D1<-D1
   rval$D2<-D2  
   }
  else    {
   rval<- dcor.dist(D1=D1,D2=D2) 
   names(rval)< "dcor" 
   }
  return(rval)
} 
################################################################################
