library(coxme)
options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

# $Id: test3.s,v 1.4 2003/08/21 21:24:40 therneau Exp $
#
# Check out a mixed sparse/non-sparse cacluation
#
# The choice of "sparse" below forces group 1 to be non-sparse and
#   group 2-4 to be sparse.  It will also force the coefficients to be
#   in 2,3,4,1 order
# The test is much more complete if the grouping variable has 4 levels instead
#   of just 2.
#
tdata1 <- data.frame(time  =c(5,4,1,1,2,2,2,2,3, 1:8), 
                     status=c(0,1,1,0,1,1,1,0,0, rep(1:0,4)),
                     x1    =c(0,1,2,0,1,1,0,1,0, rep(4:1,2)),
                     wt    =c(1,2,1,2,3,4,3,2,1, rep(1:2,4)),
                     x2    =c(1,3,5,2,3,6,4,3,1, rep(1:2,4)),
                     grp   =c(1,1,2,2,1,1,2,2,1, rep(3,4), rep(4,4)))

theta=.77
zero <- rep(0, nrow(tdata1))
fit0 <- coxme(Surv(zero, time, status) ~ x1 + x2 + (1|grp), data=tdata1,
              vfixed=theta, weight=wt, iter=0,
              sparse=c(2, .25), ties='breslow')

tfit <- coxph(Surv(time, status) ~ I(grp==2) + I(grp==3) +I(grp==4) + 
              I(grp==1) + x1 + x2,
              data=tdata1, x=T, weight=wt, iter=0, method='breslow')
dt0 <- coxph.detail(tfit)
aeq(apply(dt0$score,2,sum), fit0$u)

h0 <- apply(dt0$imat,1:2,sum) + diag(c(1/theta, 1/theta,1/theta, 1/theta,0,0))
h0[1,2:3] <- h0[2:3,1] <- 0  #sparse part
h0[2,3] <- h0[3,2] <- 0

aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(diag(gchol(h0)), diag(fit0$hmat))
hinv <- solve(h0)
hinv[1,2:3] <-hinv[2:3,1] <- 0
hinv[2,3] <- hinv[3,2] <- 0
aeq(hinv, as.matrix(fit0$var))  

fit1 <- coxme(Surv(zero, time, status) ~ x1 + x2 + (1|grp), data=tdata1,
              vfixed=theta, weight=wt, iter=1,
              sparse=c(2, .25), ties='breslow')

temp <- solve(fit0$hmat, fit0$u)
temp[1:4] <- temp[1:4] - mean(temp[1:4])
aeq(temp, c(unlist(ranef(fit1)), fixef(fit1)))

# Now test out Breslow/cox
fit0 <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata1,
              vfixed=theta, weight=wt, iter=0,
              sparse=c(2, .25), ties='breslow')
aeq(apply(dt0$score,2,sum), fit0$u)
aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(diag(gchol(h0)), diag(fit0$hmat))
aeq(hinv, as.matrix(fit0$var))  

#Efron/Cox
fit0 <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata1,
              vfixed=theta, weight=wt, iter=0,
              sparse=c(2, .25))

tfit <- coxph(Surv(time, status) ~ I(grp==2) + I(grp==3) +I(grp==4) + 
              I(grp==1) + x1 + x2,
              data=tdata1, x=T, weight=wt, iter=0)
dt0 <- coxph.detail(tfit)
aeq(apply(dt0$score,2,sum), fit0$u)

h0 <- apply(dt0$imat,1:2,sum) + diag(c(1/theta, 1/theta,1/theta, 1/theta,0,0))
h0[1,2:3] <- h0[2:3,1] <- 0  #sparse part
h0[2,3] <- h0[3,2] <- 0

aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(diag(gchol(h0)), diag(fit0$hmat))
hinv <- solve(h0)
hinv[1,2:3] <-hinv[2:3,1] <- 0
hinv[2,3] <- hinv[3,2] <- 0
aeq(hinv, as.matrix(fit0$var))  

# Efron/ag
fit0 <- coxme(Surv(zero,time, status) ~ x1 + x2 + (1|grp), data=tdata1,
              vfixed=theta, weight=wt, iter=0,
              sparse=c(2, .25))
aeq(apply(dt0$score,2,sum), fit0$u)
aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(diag(gchol(h0)), diag(fit0$hmat))
aeq(hinv, as.matrix(fit0$var))  

rm(tfit, dt0, fit0, h0, hinv, fit1, theta, zero)
