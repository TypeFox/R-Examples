set.seed(1234)

library(hgam)

test.d <- dgp(1000)
test.m <- hgam(y ~ ., data = test.d)

print(test.m)
coef(test.m)
fitted(test.m)
predict(test.m)

dgp <- function (n, sd = 1) 
{
    x1 <- seq(from = 0, to = 2 * pi, length = n)
    x2 <- runif(n = n)
    y <- rnorm(n = n, mean = sin(x1) + 2 * x2, sd = sd)
    return(data.frame(y, x1, x2))
}

sd <- 0.5
n <- 1000
dftrain <- dgp(n, sd = sd)
dftest <- dgp(10000, sd = sd)
dfplot <- dgp(100, sd = sd)

mod <- hgam(y ~ ., data = dftrain, lambda2 = 5, lambda1 = 1)
sd(dftest$y - sin(dftest$x1) - 2 * dftest$x2)
sd(dftest$y - as.numeric(predict(mod, newdata = dftest)))

layout(matrix(1:3, nr = 1))
plot(dfplot$x1, predict(mod, newdata = dfplot, which = "x1"))
lines(dfplot$x1, sin(dfplot$x1), col = "red")
plot(dfplot$x2, predict(mod, newdata = dfplot, which = "x2"))
lines(dfplot$x2, 2 * dfplot$x2, col = "red")
plot(sin(dftest$x1) + 2 * dftest$x2, predict(mod, newdata = dftest), 
     ylim = range(dftest$y))
abline(a = 0, b = 1)

n <- nrow(dftrain)
w <- rmultinom(1, n, rep(1, n) / n)
mod <- hgam(y ~ ., data = dftrain, weights = w, 
            lambda2 = 5, lambda1 = 1)

y <- mod$y
f <- fitted(mod)
w <- mod$weights
mod$model@nloglik(y, f, w) ### Risiko Trainingsdaten
mod$model@nloglik(y[w == 0], f[w == 0], rep(1, sum(w == 0))) ### Risiko Testdaten
