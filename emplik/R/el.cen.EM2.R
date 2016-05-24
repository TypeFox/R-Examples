el.cen.EM2 <- function(x,d,xc=1:length(x),fun,mu,maxit=25,error=1e-9,...){
####
#### xc is collaps control: if for index i and j, (x[] d[]) are identical
####     should they be merged into one obs. with weight 2?
####     if xc[i] != xc[j] then they do not merge. 
####     For regression applications, you should not merge -- use
####     xc = 1:length(x) .
####     The order of censoring (re-distribute weights)
####     is based on x, not fun(x).
####
   xvec <- as.vector(x) 
   d <- as.vector(d) 
   mu <- as.vector(mu)
   xc <- as.vector(xc)
   n <- length(d) 

   if (length(xvec)!=n) stop ("length of d and x must agree")
   if (length(xc)!=n) stop ("length of xc and d must agree")
   if(n <= 2*length(mu)+1) stop ("Need more observations")
   if(any((d!=0)&(d!=1)&(d!=2)))
     stop("d must be 0(right-censored) or 1(uncensored) or 2(left-censored)")
   if(!is.numeric(xvec)) stop("x must be numeric")
   if(!is.numeric(mu)) stop("mu must be numeric")

   funx <- as.matrix(fun(xvec, ...))   ##get the matrix before sorting xvec
   pp <- ncol(funx)
   if(length(mu)!=pp) stop("length of mu and ncol of fun(x) must agree")

   temp <- Wdataclean5(z=xvec,d,zc=xc, xmat=funx)  ## collaps control zc
   x <- temp$value
   d <- temp$dd
   w <- temp$weight
   funx <- temp$xxmat     ## make sure funx and x are in same order

   ###### change the last obs. among d=1,0, so that d=1
   ###### change the first obs. among d=1,2 so that d=1
   ###### this ensures we got a proper CDF for NPMLE. (no mass escapes)
   INDEX10 <- which(d != 2)
   d[ INDEX10[length(INDEX10)] ] <- 1
   INDEX12 <- which(d != 0)
   d[ INDEX12[1] ] <- 1

   xd1 <- x[d==1]
    if(length(xd1) <= 1) stop("need more distinct uncensored obs.")

   funxd1 <- funx[d==1,] 
   xd0 <- x[d==0]
   xd2 <- x[d==2]
   wd1 <- w[d==1]
   wd0 <- w[d==0]
   wd2 <- w[d==2]
   m <- length(xd0)
   mleft <- length(xd2)

##############################################
#### do the computation in 4 different cases.#
##############################################
 if( (m>0) && (mleft>0) ) { 
   pnew <- el.test.wt2(x=funxd1, wt=wd1, mu=mu)$prob
   n <- length(pnew)
   k <- rep(NA, m)
   for(i in 1:m) { k[i] <- 1+n - sum( xd1 > xd0[i] ) }
   kk <- rep(NA, mleft)
   for(j in 1:mleft) { kk[j] <- sum( xd1 < xd2[j] ) }
   num <- 1
   while(num < maxit) {
     wd1new <- wd1
     sur <- cumsumsurv(pnew)  ### rev(cumsum(rev(pnew))) 3/2015 MZ
     cdf <- 1 - c(sur[-1],0)
     for(i in 1:m)
        {wd1new[k[i]:n] <- wd1new[k[i]:n] + wd0[i]*pnew[k[i]:n]/sur[k[i]]}
     for(j in 1:mleft) {
     wd1new[1:kk[j]] <- wd1new[1:kk[j]] + wd2[j]*pnew[1:kk[j]]/cdf[kk[j]]}
     temp8 <- el.test.wt2(x=funxd1, wt=wd1new, mu=mu)
     pnew <- temp8$prob
     lam <- temp8$lambda
     num <- num +1
     }
   logel <- sum(wd1*log(pnew)) + sum(wd0*log(sur[k])) + sum(wd2*log(cdf[kk]))
   logel00 <- WCY(x=x,d=d,wt=w)$logEL   ####  should I add zc?
   }
  if( (m>0) && (mleft==0) ) {
   pnew <- el.test.wt2(x=funxd1, wt=wd1, mu=mu)$prob
   n <- length(pnew)
   k <- rep(NA,m)
   for(i in 1:m) { k[i] <- 1 + n - sum( xd1 > xd0[i] ) }
   num <- 1
   while(num < maxit) {
     wd1new <- wd1
     sur <- cumsumsurv(pnew)   ## rev(cumsum(rev(pnew))) 3/2015 MZ
     for(i in 1:m)
        {wd1new[k[i]:n] <- wd1new[k[i]:n] + wd0[i]*pnew[k[i]:n]/sur[k[i]]}
     temp9 <- el.test.wt2(x=funxd1, wt=wd1new, mu=mu)
     pnew <- temp9$prob
     lam <- temp9$lambda
     num <- num+1
     }
   sur <- cumsumsurv(pnew)   ## rev(cumsum(rev(pnew))) 3/2015 MZ
   logel <- sum( wd1*log(pnew)) + sum( wd0*log(sur[k]) )
   logel00 <- WKM(x=x,d=d,zc=xc,w=w)$logel
   }
  if( (m==0) && (mleft>0) ) {
   kk <- rep(NA, mleft)
   for(j in 1:mleft) { kk[j] <- sum( xd1 < xd2[j] ) }
   pnew <- el.test.wt2(x=funxd1, wt=wd1, mu=mu)$prob
   n <- length(pnew)
   num <- 1
   while(num < maxit) {
     wd1new <- wd1
     cdf <- cumsum(pnew) 
     for(j in 1:mleft)
       {wd1new[1:kk[j]] <- wd1new[1:kk[j]] + wd2[j]*pnew[1:kk[j]]/cdf[kk[j]]}
     temp7 <- el.test.wt2(x=funxd1, wt=wd1new, mu=mu)
     pnew <- temp7$prob
     lam <- temp7$lambda
     num <- num+1
     }
   logel <- sum( wd1*log(pnew)) + sum( wd2*log( cdf[kk] ) )
   dleft <- d 
   dleft[dleft==2] <- 0
   templeft <- WKM(x= -x, d=dleft, zc=xc, w=w)  ### bug fix 3/2008
   logel00 <- templeft$logel   ### ?? do I need a left WKM ??
   }
  if( (m==0) && (mleft==0) ) { 
    num <- 0
    temp6 <- el.test.wt2(x=funxd1, wt=wd1, mu)
    pnew <- temp6$prob
    lam <- temp6$lambda
    logel <- sum( wd1*log(pnew) ) 
    logel00 <- sum( wd1*log( wd1/sum(wd1) ) )
  }
  tval <- 2*(logel00 - logel)
  list(loglik=logel, times=xd1, prob=pnew, lam = lam, iters=num, 
             "-2LLR"=tval, Pval= 1-pchisq(tval, df=length(mu)) )
}
