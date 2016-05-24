npower.lnorm <-
function(n=NA,power=NA,fstar=1,p=0.95,gamma=0.95){
#   For random sample from lognormal distribution
#   given Ho: Xp > Lp  at the alp = 1 - gam significance level
#         where Xp is 100p percentile and Lp is specified limit
#         Ho: Fp > 100*(1-p) 
#         where Fp is the percent of Xs > Lp
#   Reject Ho if Uf(Lp,gam) < Lp OR UX(p,gam)< Lp 
#   find exact sample size n to provide power of at least (1-beta)
#   when the true value F is F* (fstar)
#  
# USAGE: npower.lnorm(n,power,fstar,
# ARGUMENTS:
#   fstar is true percent of Xs > Lp
#   power = power of test = (1 -beta)
#   p specifies 100pth percentile of X distribution
#   gam desired confidence level = 1 - alp
# VALUE: n sample size
# NOTE: Based on non-central t distribution see 
#     Lyles and Kupper (1996) JAIHA vol 57 6-15 Equation 5
# Turn off warning message that may occur in uniroot
# calls to pt and does not effect the precision of the final result
options(warn=-1) ; eps <- 0.0001

if( fstar > 100*(1-p) ) 
stop(paste("fstar must be less than",100*(1-p),"percent") )
if( is.na(n) )
 { if( is.na(fstar)||is.na(p)||is.na(gamma))stop("Missing fstar p or gamma")
   if( is.na(power) || power < eps || power > (1-eps) )
   stop("power is missing or invalid")
# local function np2 used by uniroot to find n
np2<-function(n,fs=fstar,pow=power,pv=p,ga=gamma){
ncp0<- -sqrt(n)*qnorm(pv)
ncp1<- -sqrt(n)*qnorm( 1 - fs/100 )
# t0<- tinv(n,ncp0,1-ga)
  t0<- qt(1-ga,n-1,ncp0)
#t1<- tinv(n,ncp1,pow)
  t1<- qt(pow,n-1,ncp1)
val<- (t0 - t1)
}
n <- floor(uniroot(np2,c(3,20000))$root+1  )
}
#   find POWER given true fvalue of F is f* (fstar)
# NOTE: Based on non-central t distribution see 
#     Lyles and Kupper (1996) JAIHA vol 57 6-15 Equation 4
if( is.na(power) ){
   if( is.na(fstar)||is.na(p)||is.na(gamma))stop("Missing fstar p or gamma")
   if( n < 2 || n > 20000  ) stop("n is invalid")
ncp0<- -sqrt(n)*qnorm(p)
ncp1<- -sqrt(n)*qnorm( 1 - fstar/100 )
#tvo<- tinv(n,ncp0,1-gam)
tvo<- qt(1-gamma,n-1,ncp0)
power<- pt(tvo,n-1,ncp1)
}
out<-c(n,power,fstar,p,gamma)
options(warn=0)
out
}

