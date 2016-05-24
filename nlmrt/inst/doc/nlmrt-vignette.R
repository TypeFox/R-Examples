### R code from vignette source 'nlmrt-vignette.Rnw'

###################################################
### code chunk number 1: chunk01
###################################################
options(width=60)
pastured <- data.frame(
time=c(9, 14, 21, 28, 42, 57, 63, 70, 79),
yield= c(8.93, 10.8, 18.59, 22.33, 39.35, 
         56.11, 61.73, 64.62, 67.08))
regmod <- "yield ~ t1 - t2*exp(-exp(t3+t4*log(time)))"
ones <- c(t1=1, t2=1, t3=1, t4=1) # all ones start
huetstart <- c(t1=70, t2=60, t3=0, t4=1)
require(nlmrt)


###################################################
### code chunk number 2: chunk02
###################################################
anmrt <- nlxb(regmod, start=ones, trace=FALSE, data=pastured)
print(anmrt)


###################################################
### code chunk number 3: chunk03
###################################################
anmrtx <- try(nlxb(regmod, start=huetstart, trace=FALSE, data=pastured))
print(strwrap(anmrtx))


###################################################
### code chunk number 4: chunk04
###################################################
anls <- try(nls(regmod, start=ones, trace=FALSE, data=pastured))
print(strwrap(anls))


###################################################
### code chunk number 5: chunk05
###################################################
anlsx <- try(nls(regmod, start=huetstart, trace=FALSE, data=pastured))
print(strwrap(anlsx))


###################################################
### code chunk number 6: chunk06
###################################################
awnls <- wrapnls(regmod, start=ones, data=pastured, control=list(rofftest=FALSE))
print(awnls)
cat("Note that the above is just the nls() summary result.\n")


###################################################
### code chunk number 7: chunk07
###################################################
shobbs.res <- function(x){ # scaled Hobbs weeds problem -- residual
# This variant uses looping
    if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
    y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948,
         75.995, 91.972)
    tt <- 1:12
    res <- 100.0*x[1]/(1+x[2]*10.*exp(-0.1*x[3]*tt)) - y
}
 
shobbs.jac <- function(x) { # scaled Hobbs weeds problem -- Jacobian
    jj <- matrix(0.0, 12, 3)
    tt <- 1:12
    yy <- exp(-0.1*x[3]*tt) # We don't need data for the Jacobian
    zz <- 100.0/(1+10.*x[2]*yy)
    jj[tt,1]  <-  zz
    jj[tt,2]  <-  -0.1*x[1]*zz*zz*yy
    jj[tt,3]  <-  0.01*x[1]*zz*zz*yy*x[2]*tt
    return(jj)
}


###################################################
### code chunk number 8: chunk08
###################################################
st <- c(b1=1, b2=1, b3=1)
ans1 <- nlfb(st, shobbs.res, shobbs.jac, trace=FALSE)
print(ans1)


###################################################
### code chunk number 9: chunk09
###################################################
cat("No jacobian function -- use internal approximation\n")
ans1n <- nlfb(st, shobbs.res, trace=FALSE, control=list(watch=FALSE)) # NO jacfn
print(ans1n)


###################################################
### code chunk number 10: chunk10
###################################################
shobbs.f <- function(x){
   res <- shobbs.res(x)
   as.numeric(crossprod(res))
}
shobbs.g <- function(x){
   res <- shobbs.res(x) # This is NOT efficient -- we generally have res already calculated
   JJ <- shobbs.jac(x)
   2.0*as.vector(crossprod(JJ,res))
}
require(optimx)
aopx <- optimx(st, shobbs.f, shobbs.g, control=list(all.methods=TRUE))
summary(aopx)
cat("\nNow with numerical gradient approximation or derivative free methods\n")
aopxn <- optimx(st, shobbs.f, control=list(all.methods=TRUE))
summary(aopxn) # no file output


###################################################
### code chunk number 11: chunk12
###################################################
# jres <- model2resfun(regmod, ones, funname="myxres", file="testresfn.R")
jres <- model2resfun(regmod, ones)
print(jres)
valjres <- jres(ones, yield=pastured$yield, time=pastured$time)
cat("valjres:")
print(valjres)


###################################################
### code chunk number 12: chunk13
###################################################
jjac <- model2jacfun(regmod, ones)
print(jjac)
# Note that we now need some data!
valjjac <- jjac(ones, yield=pastured$yield, time=pastured$time)
cat("valjac:")
print(valjjac)
# Now compute the numerical approximation
require(numDeriv)
Jn <- jacobian(jres, ones, , yield=pastured$yield, time=pastured$time)
cat("maxabsdiff=",max(abs(Jn-valjjac)),"\n")


###################################################
### code chunk number 13: chunk14
###################################################
ssfn <- model2ssfun(regmod, ones) # problem getting the data attached!
print(ssfn)
valss <- ssfn(ones, yield=pastured$yield, time=pastured$time)
cat("valss: ",valss,"\n")
grfn <- model2grfun(regmod, ones) # problem getting the data attached!
print(grfn)
valgr <- grfn(ones, yield=pastured$yield, time=pastured$time)
cat("valgr:")
print(valgr)
gn <- grad(ssfn, ones, yield=pastured$yield, time=pastured$time)
cat("maxabsdiff=",max(abs(gn-valgr)),"\n")


###################################################
### code chunk number 14: chunk15
###################################################
cat("\n\nHuetstart:")
print(huetstart)
valjres <- jres(huetstart, yield=pastured$yield, time=pastured$time)
cat("valjres:")
print(valjres)
valss <- ssfn(huetstart, yield=pastured$yield, time=pastured$time)
cat("valss:", valss, "\n")
valjjac <- jjac(huetstart, yield=pastured$yield, time=pastured$time)
cat("valjac:")
print(valjjac)
Jn <- jacobian(jres, huetstart, , yield=pastured$yield, time=pastured$time)
cat("maxabsdiff=",max(abs(Jn-valjjac)),"\n")
valgr <- grfn(huetstart, yield=pastured$yield, time=pastured$time)
cat("valgr:")
print(valgr)
gn <- grad(ssfn, huetstart, yield=pastured$yield, time=pastured$time)
cat("maxabsdiff=",max(abs(gn-valgr)),"\n")


###################################################
### code chunk number 15: chunk16
###################################################
cat("All ones to start\n")
anlfb <- nlfb(ones, jres, jjac, trace=FALSE, yield=pastured$yield, time=pastured$time)
print(strwrap(anlfb))
cat("Huet start\n")
anlfbh <- nlfb(huetstart, jres, jjac, trace=FALSE, yield=pastured$yield, time=pastured$time)
print(strwrap(anlfbh))


###################################################
### code chunk number 16: uweeds01
###################################################
require(nlmrt)
# Data for Hobbs problem
ydat <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972)
tdat <- 1:length(ydat)
weeddata1 <- data.frame(y=ydat, tt=tdat)
start1 <- c(b1=1, b2=1, b3=1) # name parameters for nlxb, nls, wrapnls.
eunsc <-  y ~ b1/(1+b2*exp(-b3*tt))
anlxb1 <- try(nlxb(eunsc, start=start1, data=weeddata1))
print(anlxb1)


###################################################
### code chunk number 17: bweeds01
###################################################
# WITH BOUNDS
startf1 <- c(b1=1, b2=1, b3=.1) # a feasible start when b3 <= 0.25
anlxb1 <- try(nlxb(eunsc, start=startf1, lower=c(b1=0, b2=0, b3=0), 
      upper=c(b1=500, b2=100, b3=5), data=weeddata1))
print(anlxb1)


###################################################
### code chunk number 18: bweeds02
###################################################
anlsb1 <- try(nls(eunsc, start=startf1, lower=c(b1=0, b2=0, b3=0), 
     upper=c(b1=500, b2=100, b3=5), data=weeddata1, algorithm='port'))
print(anlsb1)


###################################################
### code chunk number 19: bweeds03
###################################################
## Uncon solution has bounds ACTIVE. Infeasible start
anlxb2i <- try(nlxb(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), 
           upper=c(b1=500, b2=100, b3=.25), data=weeddata1))
print(anlxb2i)
anlsb2i <- try(nls(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), 
           upper=c(b1=500, b2=100, b3=.25), data=weeddata1, algorithm='port'))
print(anlsb2i)


###################################################
### code chunk number 20: bweeds04
###################################################
## Uncon solution has bounds ACTIVE. Feasible start
anlxb2f <- try(nlxb(eunsc, start=startf1, lower=c(b1=0, b2=0, b3=0), 
   upper=c(b1=500, b2=100, b3=.25), data=weeddata1))
print(anlxb2f)
anlsb2f <- try(nls(eunsc, start=startf1, lower=c(b1=0, b2=0, b3=0), 
   upper=c(b1=500, b2=100, b3=.25), data=weeddata1, algorithm='port'))
print(anlsb2f)


###################################################
### code chunk number 21: mweeds01
###################################################
## TEST MASKS
anlsmnqm <- try(nlxb(eunsc, start=start1, lower=c(b1=0, b2=0, b3=0), 
   upper=c(b1=500, b2=100, b3=5), masked=c("b2"), data=weeddata1))
print(anlsmnqm) # b2 masked
an1qm3 <- try(nlxb(eunsc, start=start1, data=weeddata1, masked=c("b3")))
print(an1qm3) # b3 masked 
# Note that the parameters are put in out of order to test code.
an1qm123 <- try(nlxb(eunsc, start=start1, data=weeddata1, masked=c("b2","b1","b3")))
print(an1qm123) # ALL masked - fails!!


###################################################
### code chunk number 22: bmweeds01
###################################################
## BOUNDS and MASK
an1qbm2 <- try(nlxb(eunsc, start=startf1, data=weeddata1, 
    lower=c(0,0,0), upper=c(200, 60, .3), masked=c("b2")))
print(an1qbm2)
an1qbm2x <- try(nlxb(eunsc, start=startf1, data=weeddata1, 
    lower=c(0,0,0), upper=c(48, 60, .3), masked=c("b2")))
print(an1qbm2x)


###################################################
### code chunk number 23: bmweeds10
###################################################
hobbs.res <- function(x){ # Hobbs weeds problem -- residual
    if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
    y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 50.156, 62.948,
         75.995, 91.972)
    tt <- 1:12
    res <- x[1]/(1+x[2]*exp(-x[3]*tt)) - y
}
 
hobbs.jac <- function(x) { # Hobbs weeds problem -- Jacobian
    jj <- matrix(0.0, 12, 3)
    tt <- 1:12
    yy <- exp(-x[3]*tt)
    zz <- 1.0/(1+x[2]*yy)
    jj[tt,1]  <-  zz
    jj[tt,2]  <-  -x[1]*zz*zz*yy
    jj[tt,3]  <-  x[1]*zz*zz*yy*x[2]*tt
    return(jj)
}
# Check unconstrained
ans1 <- nlfb(start1, hobbs.res, hobbs.jac)
ans1
## No jacobian - use internal approximation
ans1n <- nlfb(start1, hobbs.res) 
ans1n
# Bounds -- infeasible start
ans2i <- try(nlfb(start1, hobbs.res, hobbs.jac, 
   lower=c(b1=0, b2=0, b3=0), upper=c(b1=500, b2=100, b3=.25)))
ans2i
# Bounds -- feasible start
ans2f <- nlfb(startf1, hobbs.res, hobbs.jac, 
   lower=c(b1=0, b2=0, b3=0), upper=c(b1=500, b2=100, b3=.25))
ans2f
# Mask b2
ansm2 <- nlfb(start1, hobbs.res, hobbs.jac, maskidx=c(2))
ansm2
# Mask b3
ansm3 <- nlfb(start1, hobbs.res, hobbs.jac, maskidx=c(3))
ansm3
# Mask all -- should fail
ansma <- try(nlfb(start1, hobbs.res, hobbs.jac, maskidx=c(3,1,2)))
ansma
# Bounds and mask
ansmbm2 <- nlfb(startf1, hobbs.res, hobbs.jac, maskidx=c(2),
      lower=c(0,0,0), upper=c(200, 60, .3))
ansmbm2
# Active bound
ansmbm2x <- nlfb(startf1, hobbs.res, hobbs.jac, maskidx=c(2),
      lower=c(0,0,0), upper=c(48, 60, .3))
ansmbm2x


###################################################
### code chunk number 24: vmcgcheck
###################################################
require(Rcgmin)
require(Rvmmin)
hobbs.f <- function(x) {
   res<-hobbs.res(x)
   as.numeric(crossprod(res))
}
hobbs.g <- function(x) {
   res <- hobbs.res(x) # Probably already available
   JJ <- hobbs.jac(x)
   2.0*as.numeric(crossprod(JJ, res))
}

# Check unconstrained
a1cg <- Rcgmin(start1, hobbs.f, hobbs.g)
a1cg
a1vm <- Rvmmin(start1, hobbs.f, hobbs.g)
a1vm
## No jacobian - use internal approximation
a1cgn <- try(Rcgmin(start1, hobbs.f))
a1cgn
a1vmn <- try(Rvmmin(start1, hobbs.f))
a1vmn
# But 
grfwd <- function(par, userfn, fbase=NULL, eps=1.0e-7, ...) {
   # Forward different gradient approximation
   if (is.null(fbase)) fbase <- userfn(par, ...)  # ensure we function value at par
   df <- rep(NA, length(par))
   teps <- eps * (abs(par) + eps)
   for (i in 1:length(par)) {
      dx <- par
      dx[i] <- dx[i] + teps[i]
      df[i] <- (userfn(dx, ...) - fbase)/teps[i]
   }
   df
}
a1vmn <- try(Rvmmin(start1, hobbs.f, gr="grfwd"))
a1vmn
# Bounds -- infeasible start
# Note: These codes move start to nearest bound
a1cg2i <- Rcgmin(start1, hobbs.f, hobbs.g, 
    lower=c(b1=0, b2=0, b3=0), upper=c(b1=500, b2=100, b3=.25))
a1cg2i
a1vm2i <- Rvmmin(start1, hobbs.f, hobbs.g, 
    lower=c(b1=0, b2=0, b3=0), upper=c(b1=500, b2=100, b3=.25))
a1vm2i # Fails to get to solution!
# Bounds -- feasible start
a1cg2f <- Rcgmin(startf1, hobbs.f, hobbs.g, 
    lower=c(b1=0, b2=0, b3=0), upper=c(b1=500, b2=100, b3=.25))
a1cg2f
a1vm2f <- Rvmmin(startf1, hobbs.f, hobbs.g, 
    lower=c(b1=0, b2=0, b3=0), upper=c(b1=500, b2=100, b3=.25))
a1vm2f # Gets there, but only just!
# Mask b2
a1cgm2 <- Rcgmin(start1, hobbs.f, hobbs.g, bdmsk=c(1,0,1))
a1cgm2
a1vmm2 <- Rvmmin(start1, hobbs.f, hobbs.g, bdmsk=c(1,0,1))
a1vmm2

# Mask b3
a1cgm3 <- Rcgmin(start1, hobbs.f, hobbs.g, bdmsk=c(1,1,0))
a1cgm3
a1vmm3 <- Rvmmin(start1, hobbs.f, hobbs.g, bdmsk=c(1,1,0))
a1vmm3

# Mask all -- should fail
a1cgma <- Rcgmin(start1, hobbs.f, hobbs.g, bdmsk=c(0,0,0))
a1cgma
a1vmma <- Rvmmin(start1, hobbs.f, hobbs.g, bdmsk=c(0,0,0))
a1vmma

# Bounds and mask
ansmbm2 <- nlfb(startf1, hobbs.res, hobbs.jac, maskidx=c(2),
      lower=c(0,0,0), upper=c(200, 60, .3))
ansmbm2
a1cgbm2 <- Rcgmin(start1, hobbs.f, hobbs.g, bdmsk=c(1,0,1),
       lower=c(0,0,0), upper=c(200, 60, .3))
a1cgbm2
a1vmbm2 <- Rvmmin(start1, hobbs.f, hobbs.g, bdmsk=c(1,0,1),
       lower=c(0,0,0), upper=c(200, 60, .3))
a1vmbm2
# Active bound
a1cgm2x <- Rcgmin(start1, hobbs.f, hobbs.g, bdmsk=c(1,0,1),
       lower=c(0,0,0), upper=c(48, 60, .3))
a1cgm2x
a1vmm2x <- Rvmmin(start1, hobbs.f, hobbs.g, bdmsk=c(1,0,1),
       lower=c(0,0,0), upper=c(48, 60, .3))
a1vmm2x


###################################################
### code chunk number 25: chunk17
###################################################
require(minpack.lm)
anlslm <- nls.lm(ones, lower=rep(-1000,4), upper=rep(1000,4), jres, jjac, yield=pastured$yield, time=pastured$time)
cat("anlslm from ones\n")
print(strwrap(anlslm))
anlslmh <- nls.lm(huetstart, lower=rep(-1000,4), upper=rep(1000,4), jres, jjac, yield=pastured$yield, time=pastured$time)
cat("anlslmh from huetstart\n")
print(strwrap(anlslmh))


