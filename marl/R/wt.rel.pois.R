wt.rel.pois <-
function(x,lambda.min,lambda.max,plot=TRUE,len=100)
{
#---------------------
# Error Check
#---------------------
if(lambda.min < 0) stop("Min. Poisson mean is negative")
if(lambda.max < 0) stop("Max. Poisson mean is negative")
if(len <= 0) stop("Length of parameter must be poistive")

x <- x; lambda <- seq(lambda.min,lambda.max,len=len)
inner.fcn <- function(lambda)
{
n <- length(x);lambda.hat <- mean(x)
#**************************
wt.rel.log <- n*(lambda.hat - lambda) + log(n/lambda.hat) + n*lambda.hat*log(lambda/lambda.hat)
wt.rel <- exp(wt.rel.log)
return(wt.rel)
}

inner.fcn <- Vectorize(inner.fcn,"lambda","x")
mle.inner <- optimize(inner.fcn, lower = lambda.min, upper = lambda.max,maximum=TRUE)$maximum
auc <- integrate(inner.fcn,lambda.min,lambda.max)$value
val <- inner.fcn(lambda)

if(plot == "TRUE")plot(lambda,val,type = "l",lwd=2)
out <- list("mle"=mle.inner,"AUC" =auc,"Val" = val,"Lambda"=lambda)
return(out)
}
