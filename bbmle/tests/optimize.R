## try to reconstruct error reported by Hofert Jan Marius
## (simpler version)

Lfun <- function(x) {
  (x-5)^2
}



library(bbmle)

lb <- 6
## first try with L-BFGS-B and bounds
m1 <- mle2(Lfun,start=list(x=7),lower=6,method="L-BFGS-B")
coef(m1)
p1 <- profile(m1)
plot(p1)
(c1 <- confint(m1,quietly=TRUE))
## all  OK

m2 <- mle2(Lfun,start=list(x=7),optimizer="optimize",
           lower=lb,upper=10)
coef(m2)
p2 <- profile(m2)
(c2 <- confint(m2))
(c2 <- confint(m2))
plot(p2,show.points=TRUE)
