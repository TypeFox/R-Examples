library(coxme)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
#
# Variable shrinkage
#
ecog0 <- 1*(lung$ph.ecog==0)
ecog1 <- 1*(lung$ph.ecog==1)
ecog2 <- 1*(lung$ph.ecog==2)
ecog3 <- 1*(lung$ph.ecog==3)

fit1 <- coxph(Surv(time, status) ~ age + ridge(ecog0, ecog1, ecog2, ecog3,
                                               scale=FALSE, theta=2), lung)

fit2 <- coxme(Surv(time, status) ~ age + (ecog0+ecog1+ecog2+ecog3 |1), lung,
              vfixed=.5)

aeq(fit1$coef, c(fixef(fit2), unlist(ranef(fit2))))
indx <- c(5,1,2,3,4) #in coxme, shrinkage variables are first
all.equal(fit1$var, as.matrix(fit2$var)[indx, indx])

fit3 <- coxme(Surv(time, status) ~ age + (1|ph.ecog), lung, vfixed=.5)
all.equal(fit2$var, fit3$var)
all.equal(fit2$loglik, fit3$loglik)

fit4 <- coxme(Surv(time, status) ~ age + (1|ph.ecog), lung)
dname <- paste("ecog", 0:3, sep='')
dummy <- matrix(diag(4), 4, dimnames=list(dname,dname))
fit5 <- coxme(Surv(time, status) ~ age + (ecog0+ecog1+ecog2+ecog3 |1), lung,
              varlist=dummy)
all.equal(fit4$log, fit5$log)
aeq(fixef(fit4), fixef(fit5))
all.equal(ranef(fit4), ranef(fit5), check.attributes=FALSE) #names will differ


