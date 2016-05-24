efraction.exact <-
function(x,gam=0.95,L=NA,logx=TRUE,wpnt=FALSE) {
#  Calculate F = exceedance fraction for limit L (e.g.OEL)
#  For complete random sample size n from normal/lognormal distribution
#  
# USAGE: efraction.exact(x,gam=0.95,L,logx=TRUE) 
# ARGUMENTS:x vector of data 
#           gam = confidence level   
#           L = Limit for exceedance fraction
#           logx if logx=T use log scale
#           wpnt if true show warning 
# # NOTE: see R bug report  RP# 9171
# VALUE: estimate of F ( as percent) 
#        exact 100gam% Uf(L,gam) and Lf(L,gam)
#       for the exceedance fraction F= Pr[ x > L]
# NOTE: ( Uf(L,gam),Lf(L,gam) ) is  100*[ 1 - 2*(1-gam) ] percent  
#       Confidence Interval for F see Section 2.3
# DETAILS: R function uniroot is used to find noncentrality 
#    parameter of noncentral t distribution to calculate CLs
#    for U= (L-mu)/sigma where F= pnorm(U) see JW Eq 16 p 366
# REFERENCES: 
#      Johnson, N. L. and Welch, B. L. (1940), Applications 
#      of the Non-Central T distribution, Biometrika, 362-389
if(!wpnt)options(warn=-1)
if( is.na(L) || L <= 0 ) stop("Value of L is missing or <= 0")
n <- length(x)
#
if( logx ) {if( any(x <= 0)) stop("all data values must be positive")
            y <- log(x) ; LL<- log(L)
           }
           else { y<- x ; LL<-L }
yb <- mean(y)
sd <- sd(y)
# del is local function called by optim 
del<- function(ncp,tv=t0,df=n-1,eps=cv)
{pt(tv,df,ncp) - eps }
#
u<- (LL - yb)/sd  ; t0<- sqrt(n)*u
cv<- gam
#    use JW eq 30 to estimate delta
dap<- t0 - qnorm(gam)*( 1 + t0^2/(2*(n-1)) )^0.5
cd<- dap*c(-2,2)
while( del(cd[1],t0,n-1,cv)*del(cd[2],t0,n-1,cv)>=0 ) {cd<- cd*2 }
u2<- uniroot(del,cd )$root
#
cv<- 1 - gam
dap<- t0 - qnorm(cv)*( 1 + t0^2/(2*(n-1)) )^0.5
cd<- dap*c(-2,2)
while( del(cd[1],t0,n-1,cv)*del(cd[2],t0,n-1,cv)>=0 ) {cd<- cd*2 } 
u1<- uniroot(del,cd )$root
out<-c(u1,u2)/sqrt(n)  # Johnson and Welch eq 16
out<- 100*c( 1-pnorm(u), 1 - pnorm(out) )
if(!wpnt)options(warn=0)
out<-list("fe"=out[1],"fe.LCL"=out[2],"fe.UCL"=out[3],"L"=L,
          "gam"=gam,"Logx"= logx)
out

}

