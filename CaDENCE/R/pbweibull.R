pbweibull <-
function(q, prob, scale, shape)
{
    if(length(prob)==1) prob <- rep(prob, length(q))
    if(length(scale)==1) scale <- rep(scale, length(q))
    if(length(shape)==1) shape <- rep(shape, length(q))
    p <- 1-prob
    p[q>0] <- 1-prob[q>0]+prob[q>0]*pweibull(q[q>0], scale=scale[q>0],
                                             shape=shape[q>0])
    p 
}

