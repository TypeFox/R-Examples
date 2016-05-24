#
# Check Deming regression on a simple data set.
#
# When var(y_i) = k var(x_i) = constant there is a closed form solution,
#  see
#
closed <- function(x, y, k, w) {
    # if w is missing assume a vector of ones.
    n <- length(x)
    if (missing(w)) wt <- rep(1/n, length(x))
    else wt <- (1/w) / sum(1/w)
    xbar <- sum(wt * x)
    ybar <- sum(wt * y)
    xmat <- cbind(y-ybar, x-xbar)
    ss <- t(xmat) %*% diag(wt) %*% xmat
    temp <- (ss[1,1] - k*ss[2,2])
    beta <- (temp + sqrt(temp^2 + 4*k*ss[1,2]^2))/ (2 * ss[1,2])
    beta
}

require(deming)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

x <- 1:10
y <- 1:10 *2.3 + c(0, -1, 2, -3, 4, -5, 6, -7, 8, -5)/2

# The deming routine uses optimize, which defaults to a tolerance of
#  .Machine$double.eps^.25, less than the default for all.equal.
tol <- .Machine$double.eps^.25
true <- closed(x,y, 1)
dfit1 <- deming(y ~ x)
aeq(coef(dfit1)[2], true, tol=tol)

dfit2 <- deming(x ~ y)
aeq(1/true, coef(dfit2)[2], tol=tol)

# verify an alternate form: the Deming angle is that rotation which gives
#  a least-squares regression coef of 0
rotate <- function(theta, x, y) {
    data.frame(x = x*cos(theta) + y* sin(theta),
               y = y*cos(theta) - x* sin(theta))
}
lfit <- lm(y ~ x, data= rotate(atan(true), x, y))
aeq(coef(lfit)[2], 0)

