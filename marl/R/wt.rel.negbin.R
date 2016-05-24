wt.rel.negbin <-
function(x,mu.min,mu.max,plot = FALSE, len = 100) 
{
#---------------------
# Error Check
#---------------------
if(mu.min < 0) stop("Min. Neg. Bin. mean is negative")
if(mu.max < 0) stop("Max. Neg. Bin. mean is negative")
if(len <= 0) stop("Length of parameter must be positive")

mu <- seq(mu.min,mu.max,len=len)
x <- x;n <- length(x);r <- mean(x)^2/abs(var(x)-mean(x))
p.mle <- n*r/(n*r+sum(x))
mu.mle <- r*(1-p.mle)/p.mle

inner.fcn <- function(mu)
{
a <- mu + r; b = mu.mle + r
wt <- n*r*mu.mle/(mu.mle*b)
rel.likld <- log(wt) + n*r*log(b/a) + sum(x)*log(mu/mu.mle) + sum(x)*log(b/a)
return(exp(rel.likld))
}
inner.fcn <- Vectorize(inner.fcn,"mu","x")

mle.inner <- optimize(inner.fcn, lower = mu.min, upper = mu.max,maximum=TRUE)$maximum
auc <- integrate(inner.fcn,mu.min,mu.max)$value
val <- inner.fcn(mu)

if(plot == "TRUE")plot(mu,val,type = "l",lwd=2)
out <- list("mle"=mle.inner,"AUC" =auc,"mu"=mu,"Val" = val)
return(out)
}
