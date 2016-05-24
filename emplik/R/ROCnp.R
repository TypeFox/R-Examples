#### Modified by Mai Zhou, 9/2007, 10/2007.
#### For two samples of right censored data.
#### Both samples use nonparametric likelihood = EL. 
#### We test the Ho: ROC curve  R(t0) = b0
#### Or testing the (1-b0)th quantile of sample one equals to
#### the (1-t0)th quantile of sample two.
#### There could be several versions: either use optimize() or optim() to
#### minimize over the value of commom quantile (the cstar value in output)
#### and either use el.cen.EM2() or use emplikH.disc()  to compute EL.
#### It seems the el.cen.EM2() version is more stable.
####
#### input: t1; right censored times, sample 1.
####        d1; censoring status, d1=1 means uncensored.

ROCnp <- function(t1, d1, t2, d2, b0, t0) {

if ( length(b0) != 1) stop ("check length of b0")
if ( length(t0) != 1) stop ("check length of t0")
if ( b0 >=1 | b0 <=0) stop ("check the value of b0")
if ( t0 >=1 | t0 <=0) stop ("check the value of t0")

tempnp2 <- WKM(x=t2, d=d2)
place2 <- sum(tempnp2$surv >= t0)
c2 <- tempnp2$times[place2]

tempnp1 <- WKM(x=t1, d=d1)
place1 <- sum(tempnp1$surv > b0)
c1 <- tempnp1$times[place1]

if(c2 <= c1)  c1 <- tempnp1$times[place1 +1]  #### need to do this??

#llr <- function(const, t1, d1, t2, d2, b0, t0) {
#             npllik2 <- emplikH.disc(x=t2, d=d2, K=log(t0), 
#                   fun=function(x,theta){as.numeric(x <= theta)},
#                   theta=const)$"discrete.-2LLR"
#             npllik1 <- emplikH.disc(x=t1, d=d1, K=log(b0),
#                   fun=function(x,theta){as.numeric(x <= theta)},
#                   theta=const)$"discrete.-2LLR"
#       return( npllik2 + npllik1 )
#       }
#
###  distribution version
llr <- function(const, t1, d1, t2, d2, b0, t0) {
             npllik1 <- el.cen.EM2(x=t1, d=d1, 
                           fun= function(x,theta){as.numeric(x <= theta)},
                           mu=1-b0, theta=const)$"-2LLR"
             npllik2 <- el.cen.EM2(x=t2, d=d2, 
                           fun= function(x,theta){as.numeric(x <= theta)},
                           mu=1-t0, theta=const)$"-2LLR"
       return(npllik1+npllik2)
       }

##temp <- optim(par=(c2+c1)/2, fn=llr, method="L-BFGS-B", lower=min(c2,c1),
##           upper=max(c2,c1), t1=t1,d1=d1,t2=t2,d2=d2,b0=b0,t0=t0)

temp <- optimize(f=llr, lower=min(c2,c1), upper=max(c2,c1),
                 t1=t1, d1=d1, t2=t2, d2=d2, b0=b0, t0=t0)

##cstar <- temp$par
cstar <- temp$minimum

##val <- temp$value
val <- temp$objective

list("-2LLR"=val, cstar=cstar)
}
