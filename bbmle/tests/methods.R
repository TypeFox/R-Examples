library(bbmle)
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
LL <- function(ymax=15, xhalf=6)
    -sum(stats::dpois(y, lambda=ymax/(1+x/xhalf), log=TRUE))
options(digits=3)
mfit0 <- mle2(y~dpois(lambda=exp(interc)),
              start=list(interc=log(mean(y))),data=d)
mfit1 <- mle2(y~dpois(lambda=exp(loglambda)),
              start=list(loglambda=log(mean(y))),data=d)

coef(mfit0)
residuals(mfit0)
AIC(mfit0)
BIC(mfit0)
vcov(mfit0)
## fitted(mfit0)  ## fails, looks for default value
predict(mfit0)  ## FIXME: doesn't expand properly (need implicit lambda~1 formula??)
set.seed(1001)
simulate(mfit0)
anova(mfit0,mfit1)
summary(mfit0)
summary(mfit1)
