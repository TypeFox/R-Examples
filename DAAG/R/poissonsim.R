"poissonsim" <-
function (x = seq(0, 1, length=101), a = 2, b = -4, intcp.sd=NULL, slope.sd=NULL, seed=NULL)
{
if (!missing(seed)) set.seed(seed)
n <- length(x)
if (!missing(intcp.sd)) a <- a + rnorm(n,sd=intcp.sd)
if (!missing(slope.sd)) b <- b + rnorm(n,sd=slope.sd)
rate <- exp(a + b*x)
y <- rpois(n, lambda=rate)
data.frame(x=x, y=y)
}

