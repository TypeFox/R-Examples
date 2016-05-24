"logisticsim" <-
function (x = seq(0, 1, length=101), a = 2, b = -4, seed=NULL)
{
if (!missing(seed)) set.seed(seed)
n <- length(x)
logit <- a + b*x
p <- exp(logit)/(1+exp(logit))
y <- rbinom(n, size=1, prob=p)
data.frame(x=x, y=y)
}

