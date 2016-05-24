###################################################
### Chap14Start
###################################################
library(mistat)
library(car)


###################################################
### AvailDisRenewDis
###################################################
set.seed(123)

Ttf <- rgamma(50, 
              shape=2, 
              scale=100)

Ttr <- rgamma(50, 
              shape=2, 
              scale=1)

AvailEbd <- availDis(ttf=Ttf,  
                     ttr=Ttr, 
                     n=1000,
                     seed=123)


RenewEbd <- renewDis(ttf=Ttf, 
                     ttr=Ttr, 
                     time=1000, 
                     n=1000)

rm(AvailEbd, RenewEbd, Ttf, Ttr)


###################################################
### PlotCdfWeibull
###################################################
set.seed(123)

plot.ecdf(
  rweibull(n=100, 
           shape=1.5, 
           scale=100),
  main="")


###################################################
### PlotExponentialQQPlot
###################################################
set.seed(123)

qqPlot(
  rexp(n=100, 
       rate=1/5), 
  distribution="exp", 
  col.lines=1)


###################################################
### PlotWeibullQQPlot
###################################################
set.seed(123)

qqPlot(
  rweibull(n=100, 
           shape=2, 
           scale=2.5), 
  distribution="weibull", 
  col.lines=1, 
  shape=2, 
  scale=2.5)


###################################################
### Survreg01
###################################################
data(FAILTIME)

library(survival)

SuRe <- survreg(
  Surv(time=FAILTIME) ~ 1 , 
  dist = "exponential")

summary(SuRe)
confint(SuRe)


###################################################
### BootSurvregBeta
###################################################
library(boot)

FAILTIME[FAILTIME >= 7000] <- 7000 # Censor data at 7000

X <- data.frame(
  time= FAILTIME, 
  event=ifelse(FAILTIME < 7000, 
               yes=1, 
               no=0))

head(X, 8)


B <- boot(data=X, 
          statistic=function(x, i){
            coefficients(
              survreg(
                Surv(
                  time=x[i,1], 
                  event=x[i,2]) ~ 1 , 
                dist = "exponential"))
          }, 
          R = 100)

boot.ci(B, 
        conf=0.95, 
        type="perc")

rm(B)


###################################################
### BootSurvregBeta
###################################################
B <- boot(
  data=X, 
  statistic=function(x, i){
    coefficients(
      survreg(
        Surv(
          time=x[i,1], 
          event=x[i,2]) ~ 1 , 
        dist="weibull"))
  }, 
  R = 100)

boot.ci(B, 
        conf=0.95, 
        type="perc")

rm(B)


###################################################
### Chap14End
###################################################
rm(X, FAILTIME, SuRe)
detach(package:boot)
detach(package:mistat)
