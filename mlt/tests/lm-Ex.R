
library("mlt")

set.seed(29)
n <- 1000	
y <- rnorm(n, 2, 1.5)
d <- data.frame(y = y)

lin <- polynomial_basis(numeric_var("y", support = range(y)), coef = c(1, 1), ci = c(-Inf, 0))
m <- ctm(lin)

o <- mlt(m, data = d)

s2ml <- sqrt(var(y) * (n - 1) / n)
1 / coef(o)[2] - s2ml

-coef(o)[1] / coef(o)[2] - mean(y)

x <- runif(n, max = 2 * pi)
y <- rnorm(n, sin(x), .25)
d <- data.frame(y = y, x = x)

plot(x, y)

Bb <- Bernstein_basis(numeric_var("x", support = c(0, 2*pi)), order = 10, ui = "zero")
m <- ctm(lin, shift = Bb)

o <- mlt(m, data = d)
1 / coef(o)[2]
p <- predict(Bb,
             newdata = data.frame(x = x), coef = coef(o)[-(1:2)])
plot(x, y)
lines(sort(x), -p[order(x)] / coef(o)[2], lwd = 2, col = "red")

x <- runif(n, max = 2 * pi)
y <- rnorm(n, 2, 1.1 + sin(x) / 2)
d <- data.frame(y = y, x = x)

plot(x, y)

Bb <- Bernstein_basis(numeric_var("x", support = c(0, 2*pi)), order = 10)
m <- ctm(lin, interacting = Bb)

o <- mlt(m, data = d)

nd <- data.frame(x = sort(x))
layout(matrix(1:2, nr = 2))
tmp <- matrix(coef(o), nrow = 2)
plot(nd$x, predict(Bb, nd, coef = tmp[1,]))
plot(nd$x, predict(Bb, nd, coef = tmp[2,]))
lines(nd$x, 1.1 + sin(nd$x) / 2, col = "red")
