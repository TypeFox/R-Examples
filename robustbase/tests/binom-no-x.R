
library(robustbase)
### "intercept only" : "no x"


set.seed(101)
k <- rbinom(100, size=3, pr = 0.2)
y <- cbind(k, n.k = 3 - k)

gg <- glm(y ~ 1, family = "binomial")
(cfK <- coef(summary(gg)))

Inf. <- 1e5 # FIXME (note that much larger values *deteriorate* slightly!)
rg.Inf <- glmrob(y ~ 1, family = "binomial", tcc= Inf.)
stopifnot(all.equal(unname(cfK[1:2]),
		    unname(unlist(coef(summary(rg.Inf))[1:2])),
		    tolerance = 1e-7))# 4.09e-8

rg.0 <- glmrob(y ~ 1, family = "binomial")
summary(rg.0)
str(rg.0, digits= 6)
