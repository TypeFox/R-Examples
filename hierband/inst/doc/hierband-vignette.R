## ------------------------------------------------------------------------
library(hierband)
set.seed(123)
p <- 100
n <- 50
K <- 10
true <- ma(p, K)
x <- matrix(rnorm(n*p), n, p) %*% true$A
S <- cov(x)

## ---- fig.height=2,fig.width=4-------------------------------------------
par(mfrow=c(1,2),mar= rep(0.1, 4))
image(true$Sig,axes=F)
image(S, axes=F)

## ------------------------------------------------------------------------
library(hierband)
path <- hierband.path(S)

## ----fig.height=4,fig.width=4--------------------------------------------
par(mfrow = c(4, 5), mar = 0.1 + c(0, 0, 2, 0))
for (i in seq_along(path$lamlist))
  image(path$P[, , i], axes = F,
        main = sprintf("lam=%s", round(path$lamlist[i], 2)))

## ----fig.height=4,fig.width=4,tidy=FALSE---------------------------------
par(mfrow = c(4, 5), mar = 0.1 + c(0, 0, 2, 0))
for (i in seq_along(path$lamlist))
  image(path$P[,,i] != 0, axes = F,
        main = sprintf("lam=%s", round(path$lamlist[i], 2)))

## ------------------------------------------------------------------------
cv <- hierband.cv(path, x)
fit <- hierband(S, lam = cv$lam.best)
plot(path$lamlist, cv$m, main = "CV Frob Error", type="o",
     ylim = range(cv$m - cv$se, cv$m + cv$se), pch = 20)
lines(path$lamlist, cv$m + cv$se)
lines(path$lamlist, cv$m - cv$se)
abline(v = path$lamlist[c(cv$ibest, cv$i.1se.rule)], lty = 2)

## ------------------------------------------------------------------------
sqrt(mean((fit - true$Sig)^2))
sqrt(mean((S - true$Sig)^2))

