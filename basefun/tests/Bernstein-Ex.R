
library("basefun")

### check the approximation of a number of functions
f1 <- function(x) qnorm(pchisq(x, df = 3))
fun <- list(sin, cos, sqrt, log, f1)
dfun <- list(cos, function(x) -sin(x), function(x) -.5 * x^(-3/2), function(x) 1/x, 
             function(x) 1 / dnorm(qnorm(pchisq(x, df = 3))) * dchisq(x, df = 3))
### http://r.789695.n4.nabble.com/Derivative-of-the-probit-td2133341.html
ord <- 3:10
x <- seq(from = 0.01, to = 2*pi - 0.01, length.out = 100)
xvar <- numeric_var("x", support = range(x) + c(-.5, .5))
for (i in 1:length(fun)) {
    for (o in ord) {
        y <- fun[[i]](x)
        Bb <- Bernstein_basis(xvar, order = o)
        m <- lm(y ~ Bb(x) - 1, data = data.frame(y = y, x = x))
        R <- summary(m)$r.squared
        layout(matrix(1:2, ncol = 2))
        plot(x, fun[[i]](x), type = "l", col = "red", main = paste(deparse(fun[[i]]), o, R, sep = ":"))
        lines(x, fitted(m))
        plot(x, dfun[[i]](x), type = "l", col = "red", main = paste(deparse(fun[[i]]), o, R, sep = ":"))
        lines(x, predict(Bb, newdata = data.frame(x = x), deriv = c(x = 1), coef = coef(m)))
    }
}

### check linear extrapolation
order <- 50
xg <- seq(from = -1, to = 1, length.out = order + 1)
xvar2 <- numeric_var("x", support = c(-1.0, 1.0))
B <- Bernstein_basis(xvar2, order = order)
cf <- xg^2
x <- -150:150/100
X <- model.matrix(B, data = data.frame(x = x))
plot(x, x^2, type = "l", col = "red")
lines(x, X %*% cf)

### deriv is constant outside support
d <- model.matrix(B, data = data.frame(x = x), deriv = c(x = 1)) %*% cf
stopifnot(length(unique(round(abs(d[abs(x) > 1])), 3)) == 1)

### Legendre to Bernstein
## Example from doi: 10.1016/j.amc.2007.09.050
A <- basefun:::L2B(4L)
B <- cbind(1, c(-1, -.5, 0, .5, 1), 
              c(1, -.5, -1, -.5, 1), 
              c(-1, 2, 0, -2, 1),  
              c(1, -4, 6, -4, 1))
stopifnot(max(abs(A - B)) < .Machine$double.eps)

### Bernstein to Legendre
stopifnot(max(abs(solve(A) - basefun:::B2L(4L))) < sqrt(.Machine$double.eps))
