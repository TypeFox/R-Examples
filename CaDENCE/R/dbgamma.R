dbgamma <-
function(x, prob, scale, shape)
{
    if(length(prob)==1) prob <- rep(prob, length(x))
    if(length(scale)==1) scale <- rep(scale, length(x))
    if(length(shape)==1) shape <- rep(shape, length(x))
    d <- 1-prob
    d[x>0] <- prob[x>0]*dgamma(x[x>0], scale=scale[x>0], shape=shape[x>0])
    d
}
