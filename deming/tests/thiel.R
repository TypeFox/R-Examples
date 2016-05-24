#
# Tests of Thiel-Sen regresion
#
library(deming)

# Data from Sen, JASA 1968
# Note that his confidence interval is incorrect in the paper, since
#  his list of values for (y_i - y_j)/(t_i - t_j) has mistakes
# He also does not use the large sample formula, but rather a table from
#  Kendall for the width of the interval
sen <- data.frame(t=c(1,  2,  3,  4, 10, 12, 18),
                  y=c(9, 15, 19, 20, 45, 55, 78))
fit <- thielsen(y ~ t, data=sen, conf=.93)
all.equal(coef(fit), c("(Intercept)"=6, "t"=4))
all.equal(fit$ci[2,], c(26/7, 35/8))

# The data is strictly monotone, so the symmetric solution = usual one
fit2 <- thielsen(y ~ t, data=sen, symmetric=TRUE)
all.equal(coef(fit2), coef(fit))


# How do the symmetric and regular compare?
fcoef <- matrix(0, 7, 3)
for (i in 1:7) {
    tfit1 <- thielsen(new.lot ~ old.lot, ferritin, subset=(period==i))
    tfit2 <- thielsen(new.lot ~ old.lot, ferritin, subset=(period==i), 
                      symmetric=TRUE)
    tfit3 <- pbreg(new.lot ~ old.lot, ferritin, subset=(period==i))
    fcoef[i,] <- c(tfit1$coef[2], tfit2$coef[2], tfit3$coef[2])
}
# For this set Passing-Bablock and STS agree
# (When the slope is exactly 1 they must be identical, and these are close to 1).
all.equal(fcoef[,2], fcoef[,3])
