percentile.ml <-
function(dd,p=0.95,gam=0.95,dat=TRUE){
#  Calculate with confidence level gam that 
#  100*p percent of population lies  below the upper bound xpu
#     and above the lower bound xpl so that
#  (xpl,xpu) is an approximate 100*[ 1 - 2*(1-p) ] percent  
#  Confidence Interval for the Xp percentile of lognormal distribution
# USAGE: xpclnd(dd,p,gam,dat)
# ARGUMENTs: if dat=FALSE  ae is output matrix from mlndln(dd)
#            if dat=TRUE  dd is an n by 2 matrix or data frame
#            p is percentile and 
#            gam is confidence level
# VALUE: Xp is pth percentile of lognormal distribution
#        (xpl,xpu) METHOD 1 Confidence Interval for Xp
#        

#   for the pth percentile, i.e UTL-pg LARGE SAMPLE RESULT
#    If x is lognormal  y=log(x) is normal(mean=ym,sd=ysd)
#    yp = ym + zp*ysd is pth percentile of y
#    ypu = ym + t(gam,nt1-1 )*zp*sdyq  is upper bound
#    var(yqu) = var(yb) + zp^2*var(ysd) +2*zp*cov(yb,ysd)
#    zp is pth percentile of N(0,1) and t() is quantile of  t distn
# NOTE: The upper bound xpu is the UTL ( upper tolerance limit)
if(dat==TRUE) ae<-lnorm.ml(dd)
if(dat==FALSE) ae <- dd
mu <- ae$mu ; sig <- ae$sig ; m <-ae$m 
vmu <- ae$se.mu^2 ;  vsig <- ae$se.sigma^2
zt <- qnorm(p) ; yp<- ae$mu + zt*ae$sig 
#  cov(mu,sig) 
sdyp <- sqrt( vmu + zt^2*vsig + 2*zt*ae$cov )
xp <- exp(yp)
xpu <- exp( yp + qt(gam,m-1 )*sdyp )
xpl <- exp( yp - qt(gam,m-1 )*sdyp )
out<-c(xp,xpl,xpu)
out<-list("Xp"=xp,"Xp.LCL"=xpl,"Xp.UCL"=xpu,"p"=p,"gam"=gam)
out
}

