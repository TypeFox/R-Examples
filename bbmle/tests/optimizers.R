library(bbmle)
old_opts <- options(digits=3)
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)
suppressWarnings(fits <- lapply(c("optim","nlm","nlminb"),
               mle2,
               minuslogl=y~dpois(lambda=ymax/(1+x/xhalf)),
               start=list(ymax=15,xhalf=6),data=d,
               method="Nelder-Mead")) ## 'method' is ignored by nlm()/nlminb()

sapply(fits,coef)
sapply(fits,logLik)

(fit2 <-  mle2(y~dpois(lambda=25/(1+x/xhalf)),
              start=list(xhalf=5),data=d,
              lower=2,upper=8,
              optimizer="optimize"))

## gives error referring to 'interval' rather than 'upper'/'lower'
## (fit2 <-  mle2(y~dpois(lambda=25/(1+x/xhalf)),
##              start=list(xhalf=5),
##              optimizer="optimize"))
options(old_opts)
