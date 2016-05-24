pblnorm <-
function(q, prob, meanlog, sdlog)
{
    if(length(prob)==1) prob <- rep(prob, length(q))
    if(length(meanlog)==1) meanlog <- rep(meanlog, length(q))
    if(length(sdlog)==1) sdlog <- rep(sdlog, length(q))
    p <- 1-prob
    p[q>0] <- 1-prob[q>0]+prob[q>0]*plnorm(q[q>0], meanlog=meanlog[q>0],
                                           sdlog=sdlog[q>0])
    p 
}
