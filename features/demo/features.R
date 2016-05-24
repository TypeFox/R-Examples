##################################################
# Estimating the smooth and the derivatives of a noisy and discretely sampled function. 
n <- 200
x <- sort(runif(n))
y <- exp(-0.2 * sin(2*pi*x)) + rnorm(n, sd=0.05)

ans <- features(x, y, fits.return=TRUE, control=list(plot.it=TRUE))
ans$f
ans$f["fwiggle"]
ans$out
fget(ans)

fits <- attr(ans, "fits")
# fits$fn is the function value at locations fits$x
# fits$d1 is the first derivative at locations fits$x
# fits$d2 is the second derivative at locations fits$x

par(ask=TRUE)
plot(fits$x, fits$d1, type="l", xlab="x", ylab="First derivative")
yexact <- -0.2 * 2*pi * cos(2*pi*x) * exp(-0.2 * sin(2*pi*x)) 
lines(x, yexact, col=2)
legend(x=0, y=0.98*max(fits$d1), legend=c("Estimated", "Exact"), lty=1, col=1:2)

ans2 <- features(x, y, smoother="smooth.spline", fits.return=TRUE, control=list(plot.it=TRUE))
fits2 <- attr(ans2, "fits")

plot(fits2$x, fits2$d1, type="l", xlab="x", ylab="First derivative")
yexact <- -0.2 * 2*pi * cos(2*pi*x) * exp(-0.2 * sin(2*pi*x)) 
lines(x, yexact, col=2)
legend(x=0, y=0.98*max(fits2$d1), legend=c("Estimated", "Exact"), lty=1, col=1:2)

# Passing arguments to smoothers

ans3 <- features(x, y, fits.return=TRUE, control=list(plot.it=TRUE), bandwidth=0.2)
fits3 <- attr(ans3, "fits")

plot(fits3$x, fits3$d1, type="l", xlab="x", ylab="First derivative")
yexact <- -0.2 * 2*pi * cos(2*pi*x) * exp(-0.2 * sin(2*pi*x)) 
lines(x, yexact, col=2)
legend(x=0, y=0.98*max(fits3$d1), legend=c("Estimated", "Exact"), lty=1, col=1:2)

par(ask=FALSE)


