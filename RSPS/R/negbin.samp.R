negbin.samp <-
function(power,lambda1,k,disp,alpha = 0.05,seed = 20,numsim=1000,sig=3)
{
#------------------------------------------
# Error check
#------------------------------------------
if(alpha >= 1 || alpha <= 0) stop("Error: alpha must be between 0 and 1")
if(1 %in% k || any(k <= 0) ) stop("Error: k = 1 or k <= 0")
if( any(power >= 1) || any(power <= 0) ) stop("Error: power must be between 0 and 1")

#------------------------------------------
l.pow <- list(power=power,lambda1=lambda1,k=k,disp=disp,alpha=alpha,seed = seed,numsim=numsim)
store.pow <- expand.grid(l.pow)

inner.fcn.pow <- function(power,lambda1,k,disp,alpha,seed,numsim)
{
n.temp <- 2;Power <- 0.001
while (Power <= power)
{
n=n.temp;lambda1=lambda1;k=k;disp=disp;alpha = alpha;seed = seed;numsim=numsim
temp <- negbin.pow(n,lambda1,k,disp,alpha,seed,numsim,monitor=FALSE,sig=sig)
n.temp <- n.temp + 1
Power <- max(temp$Power)
}
out.samp <- negbin.pow(c(n.temp-1),lambda1,k,disp,alpha,seed,numsim,monitor=FALSE,sig=sig)
return(out.samp)
}
sample <- mdply(store.pow, inner.fcn.pow)
#------------------------------------------
# Organizing output
#------------------------------------------
sample.out <- cbind(sample[,-(1:dim(store.pow)[2])],Power.Expected=power)
sample.out <- sample.out[,c(8,2:5,1,6,7)]; colnames(sample.out)[c(6,7)] <- c("N.est","Power.est")
return(sample.out)
}
