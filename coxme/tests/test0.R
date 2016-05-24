library(coxme)
options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

# Really simple dataset -- covariate x1 is our old friend from the 
#  validation section at the back of the book.  (Well, as yet only
#  my online copy has the examples with weights).
#
tdata0 <- data.frame(time  =c(5,4,1,1,2,2,2,2,3), 
                     status=c(0,1,1,0,1,1,1,0,0),
                     x1    =c(0,1,2,0,1,1,0,1,0),
                     wt    =c(1,2,1,2,3,4,3,2,1),
                     x2    =c(1,3,5,2,3,6,4,3,1),
                     grp   =c(1,1,2,2,1,1,2,2,1))

# these 3 functions give the loglik, and the u/imat results for variable
#   x1
lfun <- function(beta, efron=T) {
    r <- exp(beta)
    a <- 7*r +3
    b <- 4*r +2
    temp1 <- 11*beta - (log(r^2 + 11*r +7) + 2*log(2*r +1))
    if (efron) temp2 <- (10/3)*(log(a+b) + log(a*2/3 +b) + log(a/3 +b))
    else       temp2 <- 10 * log(a+b)
    temp1 - temp2
    }

ufun <- function(beta, efron=T) {
    r <- exp(beta)
    a <- 7*r +3
    b <- 4*r +2
    xbar1 <- (2*r^2 + 11*r)/(r^2 + 11*r + 7)
    xbar2 <- 11*r/(11*r +5)
    xbar3 <- 2*r/(2*r +1)
    xbar2b <- (7*r*2/3 + 4*r)/(a*2/3 +b)
    xbar2c <- (7*r/3 + 4*r)/(a/3 + b)

    temp1 <- 11 - (xbar1 + 2*xbar3)
    if (efron) temp2 <- (10/3)*(xbar2 + xbar2b + xbar2c)
    else       temp2 <- 10*xbar2
    temp1 - temp2
    }

ifun <- function(beta, efron=T) {
    r <- exp(beta)
    a <- 7*r +3
    b <- 4*r +2
    xbar1 <- (2*r^2 + 11*r)/(r^2 + 11*r + 7)
    xbar2 <- 11*r/(11*r +5)
    xbar3 <- 2*r/(2*r +1)
    xbar2b <- (7*r*2/3 + 4*r)/(a*2/3 +b)
    xbar2c <- (7*r/3 + 4*r)/(a/3 + b)

    temp1 <- (4*r^2 + 11*r)/(r^2 + 11*r +7) - xbar1^2
    if (efron) temp2 <- (10/3)*((xbar2 - xbar2^2) + (xbar2b - xbar2b^2) +
                                (xbar2c- xbar2c^2))
    else       temp2 <- 10*(xbar2 - xbar2^2)
    temp1 + temp2 + 2*(xbar3 - xbar3^2)
    }

tfit <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=.5, weight=wt, iter=0,
              ties='breslow', sparse.calc=0)

aeq(tfit$loglik[1], lfun(0,F))
aeq(tfit$u[3], ufun(0,F))
aeq((solve(tfit$var))[3,3], ifun(0,F))

tfit1 <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=.5, weight=wt, iter=0,
              ties='breslow', sparse.calc=1)
aeq(tfit$u, tfit1$u)
all.equal(tfit$var, tfit1$var)
aeq(tfit$loglik, tfit1$loglik)

# Do the matrix form, using coxmeMlist
dmat <- diag(2)
dimnames(dmat) <- list(1:2, 1:2)
tfit2 <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=.5, weight=wt, iter=0,
              ties='breslow', varlist=dmat)
aeq(tfit$u, tfit2$u)
aeq(as.matrix(tfit$var), as.matrix(tfit2$var))
aeq(tfit$loglik, tfit2$loglik)

#Now the Efron approx
tfit <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=.5, weight=wt, iter=0,
              ties='efron')

aeq(tfit$loglik[3], lfun(0,T))
aeq(tfit$u[3], ufun(0,T))
aeq((solve(tfit$var))[3,3], ifun(0,T))


#An initial value other than 0--
tfit <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=.5, weight=wt, iter=0,
              ties='breslow', init=c(pi,0), sparse.calc=0)

aeq(tfit$loglik[3], lfun(pi,F))
aeq(tfit$u[3], ufun(pi,F))
aeq((solve(tfit$var))[3,3], ifun(pi,F))

# Check out that the old style code -- use of a separate random statement --
#  still works.  One day we will drop this test and the backwards
#  compatablility.
#
tfit1 <- coxme(Surv(time, status) ~ x1 + x2, data=tdata0,
              random= ~1|grp, vfixed=.5, weight=wt, iter=0,
              ties='breslow', init=c(pi,0), sparse.calc=1)
aeq(tfit$u, tfit1$u)
all.equal(tfit$var, tfit1$var)
aeq(tfit$loglik, tfit1$loglik)

# Efron approximation
tfit <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=.5, weight=wt, iter=0,
              ties='efron', init=c(pi,0), sparse.calc=0)
aeq(tfit$loglik[3], lfun(pi,T))
aeq(tfit$u[3], ufun(pi,T))
aeq((solve(tfit$var))[3,3], ifun(pi,T))

tfit1 <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=.5, weight=wt, iter=0,
              ties='efron', init=c(pi,0), sparse.calc=1)
aeq(tfit$u, tfit1$u)
all.equal(tfit$var, tfit1$var)
aeq(tfit$loglik, tfit1$loglik)

# Use (start, stop] style input
dummy <- rep(0, nrow(tdata0))
tfit <- coxme(Surv(dummy, time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=.5, weight=wt, iter=0,
              ties='efron', init=c(pi,0), sparse.calc=0)
aeq(tfit$loglik[3], lfun(pi,T))
aeq(tfit$u[3], ufun(pi,T))
aeq((solve(tfit$var))[3,3], ifun(pi,T))

tfit1 <- coxme(Surv(dummy, time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=.5, weight=wt, iter=0,
              ties='efron', init=c(pi,0), sparse.calc=1)
aeq(tfit$u, tfit1$u)
all.equal(tfit$var, tfit1$var)
aeq(tfit$loglik, tfit1$loglik)


#
# a copy of tdata0, but with several subjects broken into multiple pieces,
#  exercises the "add in & take out" portions of the code
#
rcnt <- c(3,2,1,1,1,2,2,1,3)  # rep count
tdata0b <- data.frame(time2  = c(3,4,5, 2,4, 1,1,2, 1,2, .5,2, 2, 1,2,3), 
		      time1  = c(0,3,4, 0,2, 0,0,0, 0,1, 0,.5, 0, 0,1,2),
		      status = c(0,0,0, 0,1, 1,0,1, 0,1, 0, 1, 0, 0,0,0),
                      x1     = rep(tdata0$x1, rcnt),
                      wt     = rep(tdata0$wt, rcnt),
                      x2     = rep(tdata0$x2, rcnt),
                      grp    = rep(tdata0$grp, rcnt))
tfit <- coxme(Surv(time1, time2, status) ~ x1 + x2 + (1|grp), data=tdata0b,
              vfixed=.5, weight=wt, iter=0,
              ties='efron', init=c(pi,0), sparse.calc=0)
aeq(tfit$loglik[3], lfun(pi,T))
aeq(tfit$u[3], ufun(pi,T))
aeq((solve(tfit$var))[3,3], ifun(pi,T))

tfit1 <- coxme(Surv(time1, time2, status) ~ x1 + x2 + (1|grp), data=tdata0b,
              vfixed=.5, weight=wt, iter=0,
              ties='efron', init=c(pi,0), sparse.calc=1)
aeq(tfit$u, tfit1$u)
all.equal(tfit$var, tfit1$var)
aeq(tfit$loglik, tfit1$loglik)

tfit <- coxme(Surv(time1, time2, status) ~ x1 + x2 + (1|grp), data=tdata0b,
              vfixed=.5, weight=wt, iter=0,
              ties='breslow', init=c(pi,0), sparse.calc=0)
aeq(tfit$loglik[3], lfun(pi,F))
aeq(tfit$u[3], ufun(pi,F))
aeq((solve(tfit$var))[3,3], ifun(pi,F))

tfit1 <- coxme(Surv(time1, time2, status) ~ x1 + x2 + (1|grp), data=tdata0b,
              vfixed=.5, weight=wt, iter=0,
              ties='breslow', init=c(pi,0), sparse.calc=1)
aeq(tfit$u, tfit1$u)
all.equal(tfit$var, tfit1$var)
aeq(tfit$loglik, tfit1$loglik)

rm (rcnt, tfit, tfit1, dummy)
