el.ltrc.EM <- function(y,x,d,fun=function(t){t},mu,maxit=30,error=1e-9) {
   yvec <- as.vector(y)
   xvec <- as.vector(x)
   nn <- length(xvec)
   if(nn <= 1) stop ("Need more observations")
   if(length(d)!=nn) stop("length of x and d must agree")
   if(any((d!=0)&(d!=1)))
     stop("d must be 0(right-censored) or 1(uncensored)")
   if(!is.numeric(xvec)) stop("x must be numeric")
   if(length(mu)!=1) stop("check the dim of constraint mu")

   yvec <- yvec[yvec > -Inf] 
   N <- length(yvec)
   if ( N == 0 ) {
    temp1 <- el.cen.EM(xvec,d,fun=fun,mu=mu)
    WILKS <- temp1$"-2LLR"
    pnew <- temp1$prob
   }

   temp <- Wdataclean2(xvec,d)
   x <- temp$value
   d <- temp$dd
   w <- temp$weight

   ###### change the last obs.'s censoring indicator, so that d=1
   ###### this ensures we got a proper df for NPMLE. (no mass escape)
   d[length(d)] <- 1

   xd1 <- x[d==1]
   funxd1 <- fun(xd1) 
   n <- length(xd1) 
    if(n <= 1) stop("need more distinct uncensored obs.")
   xd0 <- x[d==0]
   wd1 <- w[d==1]
   wd0 <- w[d==0]
   mright <- length(xd0)

   if ( mright == 0 ) {
      temp2 <- el.trun.test(yvec,xd1,fun=fun,mu=mu)
      WILKS <- temp2$"-2LLR"
      pnew <- temp2$NPMLEmu
   } 
  if( mright>0 ) {
   p0 <- LTRC(x,d,w,yvec)$survjump
   pnew <- p0
   k <- rep(NA,mright)
   for(i in 1:mright) { k[i] <- 1 + n - sum( xd1 > xd0[i] ) }
   indi <- function(u,v){ as.numeric(u > v) }
   uij <- outer(xd1,yvec,FUN="indi")
   num <- 1
   while(num <= maxit) {
     wd1new <- wd1
##########right censor 
     sur <- cumsumsurv(pnew)  ## rev(cumsum(rev(pnew)))  3/2015 MZ
     for(i in 1:mright)
        {wd1new[k[i]:n] <- wd1new[k[i]:n] + wd0[i]*pnew[k[i]:n]/sur[k[i]]}
#########left truncated
     pvec <- as.vector( pnew %*% uij )
     wd1new <- wd1new + as.vector(rowsum(t(pnew*(1-uij))/pvec, group=rep(1,N)))
######### weighted computation
     pnew <- el.test.wt(funxd1, wt=wd1new, mu=mu)$prob
     num <- num +1
     }
   sur <- cumsumsurv(pnew)  ## rev(cumsum(rev(pnew)))  3/2015 MZ
   pvec <- as.vector( pnew %*% uij )
   logel <- sum(wd1*log(pnew))+sum(wd0*log(sur[k]))-sum(log(pvec))
   sur0 <- cumsumsurv(p0)  ## rev(cumsum(rev(p0)))  3/2015 MZ
   pvec0 <- as.vector( p0 %*% uij )
   logel0 <- sum(wd1*log(p0))+sum(wd0*log(sur0[k]))-sum(log(pvec0))
   WILKS <- 2*(logel0 - logel)
  }
  list(times=xd1, prob=pnew, "-2LLR"= WILKS, Pval=1-pchisq(WILKS,df=1))
}
