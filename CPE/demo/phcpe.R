
### create a simple data set for testing
set.seed(199)
nn <- 1000
time <- rexp(nn)
status <- sample(0:1, nn, replace=TRUE)
covar <- matrix(rnorm(3*nn), ncol=3)
survd <- data.frame(time, status, covar)
names(survd) <- c("time","status","x1","x2","x3")

coxph.fit <- coxph(Surv(time,status)~x1+x2+x3,data=survd)

### Calculate CPE and CPE.SE
phcpe(coxph.fit)

### Calculate CPE only (needs much less time)
phcpe(coxph.fit, CPE.SE=FALSE)

#*** For unknown reason, 'coxph.fit' may need to be removed before running cph()***
rm(coxph.fit)

cph.fit <- cph(Surv(time, status)~x1+x2+x3, data=survd)

### Calculate CPE and CPE.SE
phcpe(cph.fit)

### Calculate CPE only (needs much less time)
phcpe(cph.fit, CPE.SE=FALSE)

