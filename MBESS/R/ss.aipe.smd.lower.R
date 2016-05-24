ss.aipe.smd.lower <- function(delta, conf.level, width, ...)
{
alpha <- 1-conf.level

# Initial starting value for n using the z distribution.
n.0 <- 2*(qnorm(1-alpha/2)/width)^2       

# Second starting value for n using the central t distribiton. 
n <- 2*((qt(1-alpha/2, 2*n.0-2))/width)^2

# Measures the discrepency between the inital and second starting values.
Difference <- abs(n-n.0)

while(Difference > 0) 
{
n.p <- n
n <- 2*((qt(1-alpha/2, 2*n-2))/width)^2
Difference <- abs(n - n.p) 
}
n <- ceiling(n)

# To ensure that the initial n is not too big (may happen with small noncentral values).
n <- n-5

# Initial estimate of noncentral value. 
# This is literally the theoretical t-value given delta and the initial estimate of sample size.
lambda.0 <- delta*sqrt(n/2)

# Initial confidence limits.
Limits.0 <- ci.smd(ncp=lambda.0, n.1=n, n.2=n, conf.level=1-alpha)

# Initial half-width for lower limit.
Diff.width.Lower.Bound <- abs(delta - Limits.0$Lower) - width

while(Diff.width.Lower.Bound > 0)
{
n <- n + 1
lambda <- delta*sqrt(n/2)
Limits <- ci.smd(ncp=lambda, n.1=n, n.2=n, conf.level=1-alpha)
Current.Lower.width <- abs(delta-Limits$Lower)
Diff.width.Lower.Bound <- Current.Lower.width - width
}
return(n)
}
