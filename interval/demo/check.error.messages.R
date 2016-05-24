library(interval)
# Check some of the calls that might be made and make sure the error messages 
# are appropriate
set.seed(1)
# Some calls that work for right censored data will not work here
test<-data.frame(aml,w=sample(c(0,1),23,replace=TRUE))
try(ictest(Surv(time, status) ~ x+strata(w), data=test))
try(ictest(Surv(time, status) ~ x+offset(w), data=test))
try(ictest(Surv(time, status) ~ x, rho=0, data=test))


