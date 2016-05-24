library(coxme)
options(na.action='na.exclude', contrasts=c('contr.treatment', 'contr.poly'))
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

#
# Test of fitting random slopes
#
# Simulation data with 9 institutions, strong age effects
#  and a random treatment effect
#      
set.seed(56)
n.subject <- seq(180, by=21, length=9) # number of subjects
slope <- sort(-.5 + rnorm(9, sd=.5))         # true treament effects

inst <- rep(1:9, n.subject)
n <- length(inst)
simdata <- data.frame(id=1:n, inst=inst,
                      trt= rep(0:1, length=n),
                      age= runif(n, 40, 70))
#risk goes up 30%/decade of age
simdata$hazard <- .8* exp(simdata$trt * rep(slope, n.subject) +
                          (simdata$age-55) * .03)

rtime <- function(hazard, censor=c(1,2)) {
    stime <- rexp(length(hazard), rate=hazard)
    ctime <- runif(length(hazard), censor[1], censor[2])
    list(time= pmin(stime, ctime), status=1*(stime <=ctime))
    }
temp <- rtime(simdata$hazard)
simdata$time <- temp$time
simdata$status <- temp$status

coef0 <- matrix(0., 2,9)
for (i in 1:9) {
    fit0 <- coxph(Surv(time, status) ~ age + trt, simdata,
                  subset=(inst==i))
    coef0[,i] <- fit0$coef
    }

# Several of these fits will differ in the last few digits on a 
#   64 bit vs 32 bit Intel processer.  The loglike is very
#   flat on top so tiny changes in the compute path lead to a small
#   change in the final solution.  Hence the "digits" argument.
fit0 <- coxph(Surv(time, status) ~ age + trt, simdata)
fit1 <- coxme(Surv(time, status) ~ age + trt + (1|inst), simdata)
print(fit1, rcoef=TRUE, digits=4)

fit2 <- coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata)
print(fit2, rcoef=TRUE, digits=3)

# And so will this one
fit3 <- coxme(Surv(time, status) ~ age + trt + (1|inst) + (trt|inst),simdata)
print(fit3, rcoef=TRUE, digits=3)

fit4 <- coxme(Surv(time, status) ~ age + trt + (1 +trt |inst), simdata)

sfit0 <- coxph(Surv(time, status) ~ age + trt + strata(inst), simdata)
sfit1 <- coxme(Surv(time, status) ~ age + trt + (trt|inst) + strata(inst),
               simdata)
print(sfit1, rcoef=TRUE, digits=4)

# Check that the start,stop code does the same
dummy <- runif(nrow(simdata), -4, -1)  #all start times before first event
fit4b <- coxme(Surv(dummy, time, status) ~ age + trt + (1 +trt |inst), simdata)
all.equal(fit4b$loglik, fit4$loglik)
all.equal(fit4b$coef, fit4$coef, tolerance=1e-7) # different order of internal
                                               # sums => tiny difference

#Comparison plot
y <- cbind(slope, NA, coef0[2,], fixef(sfit1)[2] + unlist(ranef(sfit1)),
           fixef(fit3)[2] + ranef(fit3)[[2]],
           fixef(fit4)[2] + ranef(fit4)[[1]][,2])
matplot(c(1, 1.5, 2:5), t(y), type='b', xaxt='n', xlab="Simulation", 
        ylab="Treatment coefficient", lty=1)
axis(1, 1:5, c("Sim", "Separate", "Strata", "Uncor", "Corr"))


#
# Now compute some things exactly
#
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
igchol <- function(x) {
    dd <- diag(x)
    ll <- as.matrix(x)
    ll %*% diag(dd) %*% t(ll)
    }

# For fit2
vtemp <- unlist(VarCorr(fit2))
names(vtemp) <- names(VarCorr(fit2))
fit2a <- coxme(Surv(time, status) ~ age + trt + (1|inst/trt), simdata,
               iter=0, vfixed=vtemp)
temp <- strata(simdata$inst, simdata$trt, sep='/', shortlabel=TRUE)
cfit <- coxph(Surv(time, status) ~ factor(temp) +factor(inst) +age + trt,
              simdata, iter=0, x=T)
dt2 <- coxph.detail(cfit)
u2 <- apply(dt2$score, 2, sum)
aeq(u2, fit2a$u)
imat2 <- apply(dt2$imat, 1:2, sum) + diag(c(rep(1/vtemp, c(18,9)),0,0))
aeq(imat2, as.matrix(igchol(fit2a$hmat)))

# For fit3
vtemp <- as.vector(unlist(VarCorr(fit3)))  #name not needed
fit3a <- coxme(Surv(time, status) ~ age + trt + (1|inst) + (trt|inst),
               simdata, iter=0, vfixed=as.list(vtemp))
cfit <- coxph(Surv(time, status) ~ factor(inst) * trt + age, simdata,
              iter=0, x=T)
dt3 <- coxph.detail(cfit)
u3 <- apply(dt3$score, 2, sum)
indx <- c(1:9, 12:20, 11, 10)
aeq(u3[indx], fit3a$u)
imat2 <- apply(dt3$imat, 1:2, sum)[indx,indx] + 
    diag(c(rep(1/vtemp, c(9,9)),0,0))
aeq(imat2, as.matrix(igchol(fit3a$hmat)))

fit3b <- coxme(Surv(time, status) ~ age + trt + (trt|inst) +(1|inst),
               simdata, iter=0, vfixed=as.list(rev(vtemp)))
aeq(fit3a$u, fit3b$u)
aeq(fit3b$imat, fit3b$imat)

#For sfit1
vtemp <- .0966
fit <- coxme(Surv(time, status) ~ age + trt + strata(inst) + (trt|inst),
               simdata, iter=0, vfixed=vtemp)
cfit <- coxph(Surv(time, status) ~ factor(inst):trt + trt+ age+ strata(inst),
              simdata, iter=0, x=T)
dt3 <- coxph.detail(cfit)
u3 <- apply(dt3$score, 2, sum)
indx <- c(3:11,2,1)
aeq(u3[indx], fit$u)

imat3 <- apply(dt3$imat,1:2, sum)[indx,indx] + diag(c(rep(1/vtemp,9),0,0))
aeq(imat3, as.matrix(igchol(fit$hmat)))

