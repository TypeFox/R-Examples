require(bbmle)
mle2a <- function(...)
  mle2(...)

mle2b <- function(...)
  mle2a(...)

## some data
d <- data.frame(x=0:10,y=c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8))
ym <- mean(d$y)

## some fits

(fit0 <- mle2(y~dpois(lambda=ymean),start=list(ymean=ym),data=d)) # okay
predict(fit0)
(fit0.2 <- mle2(y~dpois(lambda=ymean),start=list(ymean=ym),data=d,
                control=list(parscale=2))) # okay
predict(fit0.2)
(fit1 <- mle2a(y~dpois(lambda=ymean),start=list(ymean=ym),data=d)) # okay
(fit1.2 <- mle2a(y~dpois(lambda=ymean),start=list(ymean=ym),data=d,
                 control=list(parscale=2))) # FAILS
(fit1.3 <- mle2b(y~dpois(lambda=ymean),start=list(ymean=ym),data=d,
                 control=list(parscale=2))) # FAILS

### NOT WORKING:
if (FALSE) {
  predict(fit1)
  predict(fit1.2)
  predict(fit1.3)
}
