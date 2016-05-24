library(coxme)
options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

#
# Now do iterations on tdata0, checking it out with coxph.detail
#   Note that the coxph calls give some "singular matrix" warnings
#   which can be ignored
#
tdata0 <- data.frame(time  =c(5,4,1,1,2,2,2,2,3), 
                     status=c(0,1,1,0,1,1,1,0,0),
                     x1    =c(0,1,2,0,1,1,0,1,0),
                     wt    =c(1,2,1,2,3,4,3,2,1),
                     x2    =c(1,3,5,2,3,6,4,3,1),
                     grp   =c(1,1,2,2,1,1,2,2,1))
theta <- .53
fit0 <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=theta, weight=wt, iter=0)

tfit <- coxph(Surv(time, status) ~ I(grp==1) + I(grp==2) + x1 + x2,
	      data=tdata0, x=T, weight=wt, iter=0)
dt0 <- coxph.detail(tfit)

aeq(apply(dt0$score,2,sum), fit0$u)
h0 <- apply(dt0$imat,1:2,sum) + diag(c(1/theta, 1/theta,0,0))
aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(as.matrix(solve(fit0$var, full=F)), h0)

# Now iteration 1
fit1 <- coxme(Surv(time, status) ~ x1 + x2 +(1|grp), data=tdata0,
              vfixed=theta, weight=wt, iter=1)
aeq(fit0$u %*% fit0$var, c(unlist(ranef(fit1)), fixef(fit1)))
tfit <- coxph(Surv(time, status) ~ I(grp==1) + I(grp==2) + x1 + x2,
	      data=tdata0, x=T, weight=wt, iter=0,
	      init=c(unlist(ranef(fit1)), fixef(fit1)))
dt1 <- coxph.detail(tfit)

aeq(apply(dt1$score,2,sum)- c(unlist(ranef(fit1)), 0,0)/theta, fit1$u)
h1 <- apply(dt1$imat,1:2,sum) + diag(c(1/theta, 1/theta,0,0))
aeq(as.matrix(gchol(h1)), as.matrix(fit1$hmat))
aeq(as.matrix(solve(fit1$var, full=F)), h1)

# And iteration 2
fit2 <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=theta, weight=wt, iter=2)
aeq(solve(fit1$hmat, fit1$u), 
    c(unlist(ranef(fit2)), fixef(fit2)) - 
    c(unlist(ranef(fit1)), fixef(fit1)))
aeq(solve(fit1$hmat, fit1$u), fit1$variance %*% fit1$u)

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
fit0 <- coxme(Surv(time1, time2, status) ~ x1 + x2 + (1|grp), data=tdata0b,
              vfixed=theta, weight=wt, iter=0)

tfit <- coxph(Surv(time1, time2, status) ~ I(grp==1) + I(grp==2) + x1 + x2,
	      data=tdata0b, x=T, weight=wt, iter=0)
dt0 <- coxph.detail(tfit)

aeq(apply(dt0$score,2,sum), fit0$u)
h0 <- apply(dt0$imat,1:2,sum) + diag(c(1/theta, 1/theta,0,0))
aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(as.matrix(solve(fit0$var, full=F)), h0)

# Now iteration 1
fit1 <- coxme(Surv(time1, time2, status) ~ x1 + x2 + (1|grp), data=tdata0b,
              vfixed=theta, weight=wt, iter=1)
aeq(fit0$u %*% fit0$var, c(unlist(ranef(fit1)), fixef(fit1)))
tfit <- coxph(Surv(time1, time2, status) ~ I(grp==1) + I(grp==2) + x1 + x2,
	      data=tdata0b, x=T, weight=wt, iter=0,
	      init=c(unlist(ranef(fit1)), fixef(fit1)))
dt1 <- coxph.detail(tfit)

aeq(apply(dt1$score,2,sum)- c(unlist(ranef(fit1)), 0,0)/theta, fit1$u)
h1 <- apply(dt1$imat,1:2,sum) + diag(c(1/theta, 1/theta,0,0))
aeq(as.matrix(gchol(h1)), as.matrix(fit1$hmat))
aeq(as.matrix(solve(fit1$var, full=F)), h1)

# And iteration 2
fit2 <- coxme(Surv(time1, time2, status) ~ x1 + x2 + (1|grp), data=tdata0b,
              vfixed=theta, weight=wt, iter=2)
aeq(fit1$u %*% fit1$var, c(unlist(ranef(fit2)), fixef(fit2)) - 
                         c(unlist(ranef(fit1)), fixef(fit1)))

rm(fit0, fit1, dt0, dt1, tfit, h0, h1, theta)
rm(fit2)
