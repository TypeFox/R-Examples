
library("mlt")
set.seed(29)

n <- 100
lcf <- c("(Intercept)" = .5, "g2" = 1, "x" = 2)
sd <- .1
cf <- c("y" = 1, -lcf) / sd
cf <- cf[c("(Intercept)", "y", "g2", "x")]

d <- expand.grid(g = gl(2, 1), x = (1:n)/n)
set.seed(29)
d$y <- with(d, rnorm(nrow(d), 
           mean = lcf["(Intercept)"] + lcf["g2"] * (g == "2") + lcf["x"] * x, 
           sd = sd))

m <- ctm(polynomial_basis(numeric_var("y", support = range(d$y) * c(1, 1.1)),
                            coef = c(TRUE, TRUE), ui = diag(2), ci = c(-Inf, 0)),
           shift = ~ g + x, data = d)
mod <- mlt(m, data = d, dofit = FALSE)
coef(mod) <- cf

tfun <- function(d) with(expand.grid(d), 
    y * cf["y"] + c(cf["(Intercept)"], sum(cf[c("(Intercept)", "g2")]))[g] + cf["x"] * x)
pfun <- function(d) with(expand.grid(d), pnorm(y, 
           mean = lcf["(Intercept)"] + lcf["g2"] * (g == "2") + lcf["x"] * x,
           sd = sd))
dfun <- function(d) with(expand.grid(d), dnorm(y, 
           mean = lcf["(Intercept)"] + lcf["g2"] * (g == "2") + lcf["x"] * x,
           sd = sd))
qfun <- function(d) with(expand.grid(d), qnorm(p,
           mean = lcf["(Intercept)"] + lcf["g2"] * (g == "2") + lcf["x"] * x,
           sd = sd))

(ny <- mkgrid(mod, 10)$y)
nd <- list(y = ny, g = gl(2, 1), x = (1:10) / 10)
ndx <- expand.grid(nd[-1])
end <- ndx
end$y <- seq(from = min(d$y), to = max(d$y), length = nrow(end))

max(abs(predict(mod, newdata = end) - model.matrix(m, data = end) %*% cf))

cf2 <- cf
cf2[1:2] <- 0
max(abs(predict(mod, newdata = end, terms = "bshifting") - model.matrix(m, data = end) %*% cf2))

p1 <- predict(mod, newdata = ndx, q = ny)

p2 <- predict(mod, newdata = nd)

max(abs(c(p1) - c(p2)))

p3 <- tfun(nd)
max(abs(c(p2) - c(p3)))

p1 <- predict(mod, newdata = nd, type = "distribution")
p2 <- pfun(nd)
max(abs(c(p1) - c(p2)))

p1 <- predict(mod, newdata = nd, type = "density")
p2 <- dfun(nd)
max(abs(c(p1) - c(p2)))

p1 <- predict(mod, newdata = nd, type = "quantile", p = 1:9 / 10, K = 50000)

min(apply(matrix(p1, nrow = 9), 2, diff))

