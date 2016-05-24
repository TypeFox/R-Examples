library(coxme)
options(na.action='na.exclude')
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

#
# Similar to test1.s, but with the rats data.  This has enough groups to
#  force sparse matrix computations.
#
femrat <- rats[rats$sex=='f',] 
contr.none <- function(n,contrasts=T) {
        if(is.numeric(n) && length(n) == 1.)
                levs <- 1.:n
        else {
                levs <- n
                n <- length(n)
        }
        contr <- array(0., c(n, n), list(levs, levs))
        contr[seq(1., n^2., n + 1.)] <- 1.
	contr
	}
options(contrasts=c('contr.none', 'contr.poly'))

theta <- pi/2

# First with no sparse
fit0 <- coxme(Surv(time, status) ~ rx + (1|litter), data=femrat,
	      vfixed=theta, iter=0, sparse=c(100, .001))
tfit <- coxph(Surv(time, status) ~ factor(litter) + rx,
              data=femrat, x=T, iter=0)
dt0 <- coxph.detail(tfit)

aeq(apply(dt0$score,2,sum), fit0$u)
h0 <- apply(dt0$imat,1:2,sum) + diag(c(rep(1/theta, 50),0))
aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(diag(gchol(h0)), diag(fit0$hmat))
aeq(diag(fit0$var), diag(solve(h0)))


# then sparse
fit0 <- coxme(Surv(time, status) ~ rx + (1|litter), data=femrat,
	      vfixed=theta, iter=0, sparse=c(20, .1))

h0[1:50,1:50] <- diag(diag(h0)[1:50])
aeq(as.matrix(gchol(h0)), as.matrix(fit0$hmat))
aeq(diag(gchol(h0)), diag(fit0$hmat))
aeq(diag(fit0$var), diag(solve(h0)))

# Now iteration 1
fit1 <- coxme(Surv(time, status) ~ rx + (1|litter), data=femrat,
              vfixed=theta, iter=1, sparse=c(10, .1))
update0 <- solve(fit0$hmat, fit0$u)
update0[1:50] <- update0[1:50] - mean(update0[1:50])
aeq(update0, c(unlist(ranef(fit1)), fixef(fit1)))
tfit <- coxph(Surv(time, status) ~ factor(litter) + rx,
              data=femrat, x=T, iter=0,
              init=c(unlist(ranef(fit1)), fixef(fit1)))
dt1 <- coxph.detail(tfit)

aeq(apply(dt1$score,2,sum)- c(unlist(ranef(fit1)), 0)/theta, fit1$u)
h1 <- apply(dt1$imat,1:2,sum) + diag(c(rep(1/theta, 50),0))
h1[1:50,1:50] <- diag(diag(h1)[1:50])
aeq(as.matrix(gchol(h1)), as.matrix(fit1$hmat))
aeq(diag(gchol(h1)), diag(fit1$hmat))
aeq(diag(fit1$var), diag(solve(h1)))


# And iteration 2
fit2 <- coxme(Surv(time, status) ~ rx + (1|litter), data=femrat,
              vfixed=theta, iter=2)

update1 <- solve(fit1$hmat, fit1$u)
update1[1:50] <- update1[1:50] - mean(update1[1:50])
aeq(update1, c(unlist(ranef(fit2)), fixef(fit2)) -
    c(unlist(ranef(fit1)), fixef(fit1)))

#
# Same computation, using a specified matrix
#
fit2b <- coxme(Surv(time, status) ~ rx + (1|litter), data=femrat,
              vfixed=theta, iter=2, varlist=bdsI(seq(1, 99, 2)))
all.equal(fit2b$u, fit2$u)
all.equal(fit2b$variance, fit2$variance)
all.equal(fit2b$loglik, fit2$loglik)


