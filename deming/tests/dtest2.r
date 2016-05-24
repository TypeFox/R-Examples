#
# Weighted Deming puzzled me at first.  The Deming line with variance
#  weights does not lie in between the two weighted LS lines of y~x and x~y.
#
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)
require(deming)

#
# Create a data set where the true line is y= .5 x + 3
set.seed(12345)
true <- 1:10
x <- true + rnorm(10, mean=0, true/4)
y <- 0.5*true + rnorm(10, 3, true/4)

d1 <- deming(y ~ x, xstd=x, ystd=y)
d1b <- deming(x ~ y, ystd=x, xstd=y)
# The deming routine uses optimize, which defaults to a tolerance of
#  .Machine$double.eps^.25, less than the default for all.equal.
tol <- .Machine$double.eps^.25
aeq(coef(d1b), c(-coef(d1)[1], 1)/coef(d1)[2], tol=tol)

d2 <- deming(y ~ x, cv=TRUE)

# Given a proposed fit, this function gives the xy values on the
#  fitted line that are closest to each data point.
ufun <- function(coef, x, y, xvar, yvar) {
    wt <- 1/(yvar + coef[2]*xvar)
    u <- wt* (yvar*x + xvar* coef[2] *(y-coef[1]))
    list(x=u, y = coef[1] + coef[2]*u, wt=wt)
}

# Iterations of the proportional fit
#  It uses average weights to stabilize the iteration path
imat <- matrix(0, 20, 2)
xstd <- x; ystd <- y
imat[1,] <- coef(d1)
for (i in 2:20) {
    temp <- ufun(imat[i-1,], x, y,  xstd^2, ystd^2)
    xstd <- (temp$x + xstd)/2
    ystd <- (temp$y + ystd)/2
    ifit <- deming(y ~ x, xstd=xstd, ystd=ystd)
    imat[i,] <- coef(ifit)
}
    
# The fit with CV=true does not follow exactly this path, since for the
#  intermediate steps it uses different bounds for the optimize routine.
aeq(coef(d2), imat[10,], tol=.001)
