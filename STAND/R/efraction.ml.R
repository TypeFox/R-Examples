efraction.ml <-
function(dd,gam=0.95,L=5,dat=TRUE){
#  Calculate ML estimate of exceedance fraction F= Pr[ x > L]
#  and  "large sample" confidence limits 
#  for lognormal data with non-detects see Section 3.4
# USAGE: efclnd(dd,gam,L,dat)
# ARGUMENTS: if dat=TRUE  dd is an n by 2 matrix or data frame
#            if dat=FALSE  dd is matrix from mlndln(dd)
#            gam is one-sided confidence level(%) 
#            L is specifed limit ( e.g OEL)
# VALUE: ML estimate of exceedance fraction and 100gam% CLs
# NOTE: (Lf,Uf) is a 100*[1 - (1-gam/100)*2] Percent Confidence Interval
    if( dat==FALSE ) me <- dd
    if( dat==TRUE ) me <- lnorm.ml(dd)
#  
mu <- me$mu ; sig<- me$sigma ; LL <- log(L)
u<- ( LL - mu)/sig
pd1 <- (-1/sig) ; pd2 <- (mu - LL)/sig^2
# calculate var of u using method of statistical differentials
su <- sqrt( (pd1*me$se.mu)^2 + (pd2*me$se.sigma)^2 + 2*pd1*pd2*me$cov.musig )
m <- me$m  # number of non-detects
u1<- u + qt(gam,m-1)*su  ; u2<- u - qt(gam,m-1)*su
f1<- 1-pnorm(u1) ; f2<- 1- pnorm(u2) ; f<- 1- pnorm(u)
out<- 100*c(f,f1,f2) 


out<-list("f"=out[1],"f.LCL"=out[2],"f.UCL"=out[3],"L"=L,"gam"=gam)
out
}

