qbweibull <-
function(p, prob, scale, shape)
{
    if(length(prob)==1) prob <- rep(prob, length(p))
    if(length(scale)==1) scale <- rep(scale, length(p))
    if(length(shape)==1) shape <- rep(shape, length(p))
    q <- rep(0, length(p))
    cases <- p > (1-prob)
    q[cases] <- qweibull((prob[cases]+p[cases]-1)/prob[cases],
                         scale=scale[cases], shape=shape[cases])
    q
}

