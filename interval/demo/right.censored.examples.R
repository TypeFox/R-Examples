# Show that ictest and icfit may be applied to right censored data and show 
# when the results agree with 'survival' results and when they do not agree 
library(interval)
fit <- survfit(Surv(time, status) ~ x, data=aml)
fitIC<- icfit(Surv(time, status) ~ x, data=aml)
orig.par<-par()
par(mfrow=c(2,1))
# Notice that the Kaplan-Meier curve is undefined 
# for the Maintained group, since the largest observation is 
# censored, icfit plots this undefined region as a gray rectangle
plot(fit,main="From survfit")
plot(fitIC,main="From icfit",XLEG=100,YLEG=.8)

# The Score Statistic from ictest matches 
# Observed-Expected of the survdiff
# p-values are slightly different because the survdiff 
# uses asymptotic normality through for example Martingale methods
# while ictest uses permutation methods 
test1 <- survdiff(Surv(time, status) ~ x, data=aml)
test1
test1IC <- ictest(Surv(time, status) ~ x, data=aml)
test1IC
# Do Wilcoxon-Mann-Whitney type tests 
test2 <- survdiff(Surv(time, status) ~ x, rho=1, data=aml)
test2
test2IC <- ictest(Surv(time, status) ~ x, rho=1, data=aml)
test2IC




par(orig.par)

