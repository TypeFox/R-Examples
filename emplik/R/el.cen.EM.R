el.cen.EM <- function(x,d,fun=function(t){t},mu,maxit=25,error=1e-9, ...) {
   xvec <- as.vector(x)
   nn <- length(xvec)
   if(nn <= 1) stop ("Need more observations")
   if(length(d)!=nn) stop("length of x and d must agree")
   if(any((d!=0)&(d!=1)&(d!=2)))
     stop("d must be 0(right-censored) or 1(uncensored) or 2(left-censored)")
   if(!is.numeric(xvec)) stop("x must be numeric")
   if(length(mu)!=1) stop("check the dim of constraint mu")

   temp <- Wdataclean2(xvec,d)
   x <- temp$value
   d <- temp$dd
   w <- temp$weight

   ###### change the last obs. among d=1,0, so that d=1
   ###### change the first obs. among d=1,2 so that d=1
   ###### this ensures we got a proper dist. for NPMLE. (no mass escape)
   INDEX10 <- which(d != 2)
   d[ INDEX10[length(INDEX10)] ] <- 1
   INDEX12 <- which(d != 0)
   d[ INDEX12[1] ] <- 1

   xd1 <- x[d==1]
    if(length(xd1) <= 1) stop("need more distinct uncensored obs.")
   funxd1 <- fun(xd1, ...) 
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
   pnew <- el.test.wt(funxd1, wt=wd1, mu)$prob
   n <- length(pnew)
   k <- rep(NA, m)
   for(i in 1:m) { k[i] <- 1+n - sum( xd1 > xd0[i] ) }
   kk <- rep(NA, mleft)
   for(j in 1:mleft) { kk[j] <- sum( xd1 < xd2[j] ) }
   num <- 1
   while(num < maxit) {
     wd1new <- wd1
     sur <- cumsumsurv(pnew)
     cdf <- 1 - c(sur[-1],0)
     for(i in 1:m)
        {wd1new[k[i]:n] <- wd1new[k[i]:n] + wd0[i]*pnew[k[i]:n]/sur[k[i]]}
     for(j in 1:mleft)
        {wd1new[1:kk[j]] <- wd1new[1:kk[j]] + wd2[j]*pnew[1:kk[j]]/cdf[kk[j]]}
     temp8 <- el.test.wt(funxd1, wt=wd1new, mu)
     pnew <- temp8$prob
     lam <- temp8$lam
     num <- num +1
     }
   logel <- sum(wd1*log(pnew)) + sum(wd0*log(sur[k])) + sum(wd2*log(cdf[kk]))
   tempDB <- WCY(x=x,d=d,wt=w)
   logel00 <- tempDB$logEL
   funNPMLE <- sum( fun(tempDB$time)* tempDB$jump )
  }
  if( (m>0) && (mleft==0) ) {
   temp3 <- WKM(x=x,d=d,w=w)
   logel00 <- temp3$logel
   funNPMLE <- sum(funxd1 * temp3$jump[temp3$jump > 0])
#  now the iteration 
   pnew <- el.test.wt(x=funxd1, wt=wd1, mu=mu)$prob
   n <- length(pnew)
   k <- rep(NA,m)
   for(i in 1:m) { k[i] <- 1 + n - sum( xd1 > xd0[i] ) }
   num <- 1
   while(num < maxit) {
     wd1new <- wd1
     sur <- cumsumsurv(pnew)  ## rev(cumsum(rev(pnew))) 3/2015 MZ
     for(i in 1:m)
        {wd1new[k[i]:n] <- wd1new[k[i]:n] + wd0[i]*pnew[k[i]:n]/sur[k[i]]}
     temp9 <- el.test.wt(funxd1, wt=wd1new, mu)
     pnew <- temp9$prob
     lam <- temp9$lam
     num <- num +1
     }
   sur <- cumsumsurv(pnew)   ## rev(cumsum(rev(pnew)))  3/2015  MZ
   logel <- sum( wd1*log(pnew)) + sum( wd0*log(sur[k]) )
  }
  if( (m==0) && (mleft>0) ) {
   kk <- rep(NA, mleft)
   for(j in 1:mleft) { kk[j] <- sum( xd1 < xd2[j] ) }
   pnew <- el.test.wt(funxd1, wt=wd1, mu)$prob
   n <- length(pnew)
   num <- 1
   while(num < maxit) {
     wd1new <- wd1
     cdf <- cumsum(pnew) 
     for(j in 1:mleft)
        {wd1new[1:kk[j]] <- wd1new[1:kk[j]] + wd2[j]*pnew[1:kk[j]]/cdf[kk[j]]}
     temp7 <- el.test.wt(funxd1, wt=wd1new, mu)
     pnew <- temp7$prob
     lam <- temp7$lam
     num <- num +1
     }
   logel <- sum( wd1*log(pnew)) + sum( wd2*log( cdf[kk] ) )
   dleft <- d
   dleft[dleft==2] <- 0 
   templeft <- WKM(x= - x, d=dleft, w=w)  ### bug fix 3/2008
   logel00 <- templeft$logel    ### substitute a left WKM() ???
   funNPMLE <- NA 
  }
  if( (m==0) && (mleft==0) ) {
    funNPMLE <- sum( funxd1 * wd1/sum(wd1) )
    logel00 <- sum( wd1*log( wd1/sum(wd1) ) )
    temp6 <- el.test.wt(funxd1, wt=wd1, mu)
    pnew <- temp6$prob
    lam <- temp6$lam
    logel <- sum( wd1*log(pnew) ) 
  }
# get ready for exit
  tval <- 2*(logel00 - logel)
  list(loglik=logel, times=xd1, prob=pnew, funMLE=funNPMLE, lam=lam,
             "-2LLR"=tval, Pval= 1-pchisq(tval, df=1) )
}
