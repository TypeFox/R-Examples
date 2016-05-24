
library("mlt")
set.seed(29)

n <- 100
p <- 2
x <- matrix(runif(n * p), nrow = n)
beta <- c(1, -1)
y <- factor(rbinom(n, size = 1, prob = plogis(x %*% beta)))
df <- data.frame(y = y, x)

m1 <- glm(y ~ X1 + X2, data = df, family = binomial())
coef(m1)

m <- ctm(~ y, shift = ~ X1 + X2, todist = "Logis", data = df)
m2 <- mlt(m, data = df, fixed = c("y1" = 0))
coef(m2)

max(abs(coef(m1) + coef(m2)[-2]))

logLik(m1)
logLik(m2)

library("nnet")

m1 <- multinom(Species ~ ., data = iris)
coef(m1)

oiris <- iris
oiris$Species <- ordered(oiris$Species)

r <- as.basis(oiris$Species)
#r <- as.basis(~ Species, data = oiris, remove_intercept = TRUE,
#              contrasts.arg = list(Species = function(n)
#                  contr.treatment(n, base = 3)),
#              ui = diff(diag(2)), ci = rep(0, 1))

m <- ctm(r, interacting = as.basis(~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data = oiris),
           todistr = "Logis")
m2 <- mlt(m, data = oiris)
coef(m2)

s <- sort(unique(oiris$Species))

pp2 <- predict(m2, newdata = oiris, q = s, type = "density")

pp1 <- predict(m1, newdata = iris, type = "prob")

max(abs(pp1 - t(pp2)))

