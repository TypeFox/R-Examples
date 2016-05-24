library(bbmle)
library(testthat)
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
LL <- function(ymax=15, xhalf=6)
    -sum(stats::dpois(y, lambda=ymax/(1+x/xhalf), log=TRUE))
mfit0 <- mle2(y~dpois(lambda=exp(interc)),
              start=list(interc=log(mean(y))),data=d)

mfit1 <- mle2(y~dpois(lambda=exp(loglambda)),
              start=list(loglambda=log(mean(y))),data=d)
              
gfit0 <- glm(y~1,family=poisson)
expect_equal(unname(coef(mfit0)),unname(coef(gfit0)))
expect_equal(logLik(mfit0),logLik(gfit0))
expect_equal(predict(mfit0),  ## only one value for now
             unique(predict(gfit0,type="response")))

## FIXME: residuals are backwards
expect_equal(residuals(mfit0,type="response"),unname(residuals(gfit0,type="response")))
## FIXME: residuals are backwards
expect_equal(residuals(mfit0,type="pearson"),unname(residuals(gfit0,type="pearson")))


