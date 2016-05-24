poiss.samp <-
function(power,lambda1,k,alpha = 0.05,seed = 20,numsim=500,sig=3)
{
  
#------------------------------------------
# Error check
#------------------------------------------
if(alpha >= 1 || alpha <= 0) stop("Error: alpha must be between 0 and 1")
if(1 %in% k || any(k <= 0) ) stop("Error: k = 1 or k <= 0")
if( any(power >= 1) || any(power <= 0) ) stop("Error: power must be between 0 and 1")

#------------------------------------------
l.pow <- list(power=power,lambda1=lambda1,k=k,alpha=alpha,seed = seed,numsim=numsim)
store.pow <- expand.grid(l.pow)

inner.fcn.pow <- function(power,lambda1,k,alpha = 0.05,seed = 20,numsim=500)
{
n.temp <- 2;Power <- 0.001
while (Power <= power)
{
n=n.temp;lambda1=lambda1;k=k;alpha = alpha;seed = seed;numsim=numsim
temp <- poiss.pow(n,lambda1,k,alpha,seed,numsim,monitor=FALSE,sig=sig)
n.temp <- n.temp + 1
Power <- max(temp$Power)
}
out.samp <- poiss.pow(c(n.temp-1),lambda1,k,alpha,seed,numsim,monitor=FALSE,sig=sig)
return(out.samp)
}
sample <- mdply(store.pow, inner.fcn.pow)
#------------------------------------------
# Organizing output
#------------------------------------------
sample.out <- cbind(sample[,-(1:dim(store.pow)[2])],Power.Expected=power)
sample.out <- sample.out[,c(7,2:4,1,5,6)]; colnames(sample.out)[c(5,6)] <- c("N.est","Power.est")
return(sample.out)
  
}
