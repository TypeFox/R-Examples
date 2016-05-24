###############################################
############# emplikH2() ######################
###############################################

emplikH2.test <- function(x, d, y= -Inf, K, fun, 
	                tola = .Machine$double.eps^.5,...)
{
if(!is.numeric(x)) stop("x must be numeric values -- observed times")
n <- length(x) 
if( n <= 2 ) stop("Need more observations than two")
if( length(d) != n ) stop("length of x and d must agree") 
if(any((d!=0)&(d!=1))) stop("d must be 0/1's for censor/not-censor")

#temp <- summary(survfit(Surv(x,d),se.fit=F,type="fleming",conf.type="none"))
#
newdata <- Wdataclean2(x,d)
temp <- DnR(newdata$value, newdata$dd, newdata$weight, y=y)


Dtime <- temp$times         # only uncensored time?  Yes. 
risk <- temp$n.risk 
jump <- (temp$n.event)/risk

funtime <- fun(Dtime,...)
funh <- sqrt(n) * funtime/risk                      # that is Zi/sqrt(n)  
funtimeTjump <- funtime * jump 

if(jump[length(jump)] >= 1) funh[length(jump)] <- 0  #for inthaz and weights

inthaz <- function(x, ftj, fh, Konst){ sum(ftj/(1 + x * fh)) - Konst } 

diff <- inthaz(0, funtimeTjump, funh, K)

if( diff == 0 ) { lam <- 0 } else {
    step <- 0.2/sqrt(n) 
    ### if(abs(diff) > 99*log(n)*step )        ##why 99*log(n)? no reason, you 
    ### stop("given theta value is too far away from theta0") # need something.

    mini<-0
    maxi<-0 
    if(diff > 0) { 
    maxi <- step 
    while(inthaz(maxi, funtimeTjump, funh, K) > 0  && maxi < 99*log(n)*step) 
    maxi <- maxi+step 
    } 
    else { 
    mini <- -step 
    while(inthaz(mini, funtimeTjump, funh, K) < 0 && mini > - 99*log(n)*step) 
    mini <- mini - step 
    } 

    if(inthaz(mini,funtimeTjump,funh,K)*inthaz(maxi,funtimeTjump,funh,K)>0)
    stop("given theta is too far away from theta0")

    temp2 <- uniroot(inthaz,c(mini,maxi), tol = tola, 
                  ftj=funtimeTjump, fh=funh, Konst=K)  
    lam <- temp2$root 
} 

onepluslamh<- 1 + lam * funh   # this is 1 + lam Zi in Ref. 

weights <- jump/onepluslamh  #need to change last jump to 1? NO. see above

loglik <- 2*(sum(log(onepluslamh)) - sum((onepluslamh-1)/onepluslamh) )
#?is that right? YES see (3.2) in Ref. above. This is ALR, or Poisson LR.

#last <- length(jump)      #to compute loglik2, we need to drop last jump
#if (jump[last] == 1) {
#                     risk1 <- risk[-last]
#                     jump1 <- jump[-last]
#                     weights1 <- weights[-last]
#                     } else {
#                            risk1 <- risk
#                            jump1 <- jump
#                            weights1 <- weights
#                            }
#
#loglik2 <- 2*( sum(log(onepluslamh)) + 
#          sum( (risk1 -1)*log((1-jump1)/(1- weights1) ) )  ) 
# this version of LR seems to give negative value sometime???

list( "-2LLR"=loglik,  ### drop this output "-2logemLRv2"=loglik2, 
      lambda=lam/sqrt(n), times=Dtime, wts=weights, 
      nits=temp2$nf, message=temp2$message )
}
# what should be the fun() and K if I want to perform a (1-sample) 
# log-rank test?
# fun3 <- function(t1, z1) { psum( t( outer(z1, t1, FUN=">=") ) ) } 
# this is similar to the function in LogRank2() function. Need psum/2/3.
# And K = int R(t) dH(t)  = sum( H(z1) ) For example if H() is 
# exponential(0.3) then H(t) = 0.3*t, i.e. K <- sum(0.3* z1) 
# so finally a call may look like
#
# Assume z1 and d1 are the inputs:
# emlik2(z1, d1, sum(0.3* z1), 
#   fun3 <- function(t1,z){psum(t(outer(z,t1,FUN=">=") ) ) }, z=z1)
#
# Now use z1<-c(1,2,3,4,5) and d1<-c(1,1,0,1,1) we get
# emplik2(z1, d1, sum(0.25* z1),
#   fun3 <- function(t1,z){psum(t(outer(z,t1,FUN=">=") ) ) }, z=z1)
#
# with outputs that include this (and more)
# $ "-2logemLR":
# [1] 0.02204689
# This tests if the (censored) obs. z1 is from exp(0.25)

