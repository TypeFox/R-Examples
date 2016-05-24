percentile.exact <-
function(x,p=0.95,gam=0.95,logx=TRUE,wpnt=FALSE){
#    Find point estimates and confidence limits for the Xp percentile
#    complete sample from a normal/lognormal distribution
# USAGE: percentile.exact(x,p,gam,L.logx)
# ARGUMENTS:
#     x is a vector of positive lognormal data
#     p is proportion that should fall below L
#     gam is the one-sided confidence level 
#     L is the specified limit of interest
# VALUE: data.frame containing point estimates 
#     of Xp and F-L in column 1
#     100gam% lower confidence limits in column 2
#     100gam% upper confidence limits in column 3 
#      Xp    LX(p,gam)    UX(p,gam)
#      F-L   LCf(L,gam)   UCf(L,gam)
# NOTE: The combined confidence limits are an approximate
#       100*[2*gam -1] Percent Confidence Interval
#
# local function kf calcculate the K factor from Non-Central T distribution
kf <- function(n,p,gam){ -qt(1-gam,n-1,qnorm(1-p)*sqrt(n))/sqrt(n) }
kf2 <- function(n,p,gam){ qt(gam,n-1,qnorm(p)*sqrt(n))/sqrt(n) }
if(!wpnt)options(warn=-1)
n <- length(x)
#
#     The next line is an example of "error checking"
if( logx ) {if( any(x <= 0)) stop("all data values must be positive")
            y <- log(x)
           }
           else { y<- x }
# calculate mean(yb) and standard deviation(sy) 
yb <- mean(y)
sy <- sd(y)
n <- length(y)
if(p > 0.5) {ku <- kf(n,p,gam); kl <- kf(n,p,1-gam) }
  else { ku <- kf2(n,p,gam) ; kl <- kf2(n,p,1-gam) }

yp <- yb + qnorm(p)*sy   # point estimate 
lcl <-  yb + kl*sy  # exact LCL 
ucl <-  yb + ku*sy    # exact UCL
if(logx){
yp<- exp(yp)
lcl<- exp( lcl )
ucl<- exp(ucl)
}
out<-list("Xp"=yp,"Xpe.LCL"=lcl,"Xpe.UCL"=ucl,"p"=p,
      "gam"=gam,"Logx"= logx, "n"=n, "Ku"=ku,"Kl"=kl)
if(!wpnt)options(warn=0)
out
}

