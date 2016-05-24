library(gnm)
set.seed(1)

time <- c(21, 18, 33, 17, 35, 23, 43)
age <- unlist(sapply(time, seq, from = min(time)))
lowerMax <- min(age) - 1 #16
upperMin <- max(age) + 1 #44
leftSlope <- c(0.1, 0.2)
leftAdjust <- log(lowerMax - 14)
f <- as.factor(rep(1:2, each = 39))

family <- binomial(link = "cloglog")

y <- leftSlope[f] * log(age + exp(leftAdjust) - lowerMax)
y <- family$linkinv(y)

#don't test as Log not exported (N.B. fails if use gnm:::Log...)
#set.seed(1)
#test <- gnm(y ~ -1 + Mult(f, Log(offset(age - lowerMax) + Exp(1))),
#            family = binomial(link = "cloglog"))

#set.seed(1)
#test <- gnm(y ~ 0 + Nonlin(LogExcess(age, side = "left",
#                                     slopeFormula = ~ 0 + f)),
#            family = binomial(link = "cloglog"))
#coef(test)
