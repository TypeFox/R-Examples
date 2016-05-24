
## Examples from boot:
library(boot)

## Stratified resampling
diff.means <- function(d, f)
{    n <- nrow(d)
    gp1 <- 1:table(as.numeric(d$series))[1]
    m1 <- sum(d[gp1,1] * f[gp1])/sum(f[gp1])
    m2 <- sum(d[-gp1,1] * f[-gp1])/sum(f[-gp1])
    ss1 <- sum(d[gp1,1]^2 * f[gp1]) - (m1 *  m1 * sum(f[gp1]))
    ss2 <- sum(d[-gp1,1]^2 * f[-gp1]) - (m2 *  m2 * sum(f[-gp1]))
    c(m1 - m2, (ss1 + ss2)/(sum(f) - 2))
}
grav1 <- gravity[as.numeric(gravity[,2]) >= 7,]
boot(grav1, diff.means, R = 999, stype = "f", strata = grav1[,2])

## Nuclear data
nuke <- nuclear[, c(1, 2, 5, 7, 8, 10, 11)]
nuke.lm <- glm(log(cost) ~ date+log(cap)+ne+ct+log(cum.n)+pt, data = nuke)
nuke.diag <- glm.diag(nuke.lm)
nuke.res <- nuke.diag$res * nuke.diag$sd
nuke.res <- nuke.res - mean(nuke.res)
    
nuke.data <- data.frame(nuke, resid = nuke.res, fit = fitted(nuke.lm))
    
new.data <- data.frame(cost = 1, date = 73.00, cap = 886, ne = 0,
                           ct = 0, cum.n = 11, pt = 1)
new.fit <- predict(nuke.lm, new.data)
    
nuke.fun <- function(dat, inds, i.pred, fit.pred, x.pred)
{
    lm.b <- glm(fit+resid[inds] ~ date+log(cap)+ne+ct+log(cum.n)+pt,
                data = dat)
    pred.b <- predict(lm.b, x.pred)
    c(coef(lm.b), pred.b - (fit.pred + dat$resid[i.pred]))
}
    
nuke.boot <- boot(nuke.data, nuke.fun, R = 999, m = 1, 
                  fit.pred = new.fit, x.pred = new.data)
mean(nuke.boot$t[, 8]^2)
new.fit - sort(nuke.boot$t[, 8])[c(975, 25)]

nuke.boot <- boot(nuke.data, nuke.fun, R = 999, m = 1, 
                  fit.pred = new.fit, x.pred = new.data)

## A linear model fit:
n <- 500000
p <- 10
X <- matrix(rnorm(n * p), n, p)
y <- rnorm(n) + X[,1]
fit <- lm(y ~ X)
