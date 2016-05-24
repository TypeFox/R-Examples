acw <- function(xxx) {
xbar <- mean(xxx$x)
ybar <- mean(xxx$y)
theta   <- atan2(xxx$y - ybar,xxx$x-xbar)
theta   <- ifelse(theta > 0, theta, theta + 2 * pi)
theta.0 <- sort(unique(theta))
iii     <- match(theta.0, theta)
xxx$x   <- xxx$x[iii]
xxx$y   <- xxx$y[iii]
xxx$bp  <- xxx$bp[iii]
xxx
}
