ss.aipe.smd.full <- function(delta, conf.level, width, ...)
{
#if(warn==FALSE) options(warn=-10)
alpha <- 1-conf.level

# Initial starting value for n using the z distribution.
n.0 <- 2*(qnorm(1-alpha/2)/(width/2))^2       

# Second starting value for n using the central t distribiton. 
n <- 2*((qt(1-alpha/2, 2*n.0-2))/(width/2))^2

# measures the discrepency between the inital and second starting values.
Difference <- abs(n-n.0)

while(Difference > .000001) 
{  
n.p <- n
n <- 2*((qt(1-alpha/2, 2*n-2))/(width/2))^2
Difference <- abs(n - n.p) 
}
n <- ceiling(n)

# To ensure that the initial n is not too big.
n <- max(4, n-5)

# Initial estimate of noncentral value. 
# This is literally the theoretical t-value given delta and the initial estimate of sample size.
lambda.0 <- delta*sqrt(n/2)

# Initial confidence limits.
Limits.0 <- ci.smd(ncp=lambda.0, n.1=n, n.2=n, conf.level=1-alpha)

# Initial full-width for confidence interval.
Diff.width.Full <- abs(Limits.0$Upper - Limits.0$Lower) - width

while(Diff.width.Full > 0)
{
n <- n + 1
lambda <- delta*sqrt(n/2)
Limits <- ci.smd(ncp=lambda, n.1=n, n.2=n, conf.level=1-alpha)
Current.width <- abs(Limits$Upper - Limits$Lower)
Diff.width.Full <- Current.width - width
}

#if(warn==FALSE) options(warn=1)

return(n)
}
