library(coxme)
aeq <- function(x, y) all.equal(as.vector(x), as.vector(y))
options(na.action=na.exclude)

#
# trivial test of strata
#
tdata0 <- data.frame(time  =c(5,4,1,1,2,2,2,2,3), 
                     status=c(0,1,1,0,1,1,1,0,0),
                     x1    =c(0,1,2,0,1,1,0,1,0),
                     wt    =c(1,2,1,2,3,4,3,2,1),
                     x2    =c(1,3,5,2,3,6,4,3,1),
                     grp   =c(1,1,2,2,1,1,2,2,1))

rcnt <- c(3,2,1,1,1,2,2,1,3)  # rep count
tdata0b <- data.frame(time2  = c(3,4,5, 2,4, 1,1,2, 1,2, .5,2, 2, 1,2,3), 
		      time1  = c(0,3,4, 0,2, 0,0,0, 0,1, 0,.5, 0, 0,1,2),
		      status = c(0,0,0, 0,1, 1,0,1, 0,1, 0, 1, 0, 0,0,0),
                      x1     = rep(tdata0$x1, rcnt),
                      wt     = rep(tdata0$wt, rcnt),
                      x2     = rep(tdata0$x2, rcnt),
                      grp    = rep(tdata0$grp, rcnt))

# Undo a gchol: given gchol(x), returns x
igchol <- function(x) {
    dd <- diag(x)
    ll <- as.matrix(x)
    ll %*% diag(dd) %*% t(ll)
    }

tfit <- coxme(Surv(time, status) ~ x1 + x2 + (1|grp), data=tdata0,
              vfixed=.5, weight=wt, iter=0,
              ties='breslow', init=c(pi,0), sparse.calc=0)

sdata0 <- rbind(tdata0, tdata0)
sdata0$strat <- rep(1:2, each=nrow(tdata0))

sfit0 <- coxme(Surv(time,status) ~ x1 + x2 + strata(strat) + (1|grp), 
               data=sdata0, ties='breslow', iter=0, sparse.calc=0,
	       vfixed=0.5, weight=wt, init=c(pi,0))

sfit1 <- coxme(Surv(time,status) ~ x1 + x2 + strata(strat) + (1|grp), 
               data=sdata0, ties='breslow', iter=0, sparse.calc=1,
	       vfixed=0.5, weight=wt, init=c(pi,0))

aeq(sfit0$u, sfit1$u)
all.equal(sfit0$hmat, sfit1$hmat)
aeq(sfit0$u/2, tfit$u)
tpen <- diag(c(2,2,0,0))
aeq(igchol(sfit0$hmat)-tpen, 2*(igchol(tfit$hmat)-tpen))

#
# Repeat for start/stop data
#   The subset makes the test more rigorous (someone exits a risk
#   set without entering it immediatly again).
tfit <- coxme(Surv(time1, time2, status) ~ x1 + x2 + (1|grp), data=tdata0b,
              vfixed=0.6, weight=wt, iter=0,
              ties='efron', init=c(pi,0), subset=(-1))

sdata0b <- rbind(tdata0b[-1,], tdata0b[-1,])
sdata0b$strat <- rep(1:2, each=nrow(tdata0b)-1)

sfit0 <- coxme(Surv(time1,time2,status) ~ x1 + x2 + strata(strat) +(1|grp),
	       data=sdata0b, iter=0, sparse.calc=0,
	       vfixed=0.6, weight=wt, init=c(pi,0))
sfit1 <- coxme(Surv(time1,time2,status) ~ x1 + x2 + strata(strat) + (1|grp),
	       data=sdata0b, iter=0, sparse.calc=1,
	       vfixed=0.6, weight=wt, init=c(pi,0))
aeq(tfit$u, sfit0$u/2)
aeq(sfit0$u, sfit1$u)
tpen <- diag(c(1/.6, 1/.6, 0,0))
aeq(sfit0$imat, sfit1$imat)
aeq(igchol(sfit0$hmat)-tpen, 2*(igchol(tfit$hmat)-tpen))
