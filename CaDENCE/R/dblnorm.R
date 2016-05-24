dblnorm <-
function(x, prob, meanlog, sdlog)
{
    if(length(prob)==1) prob <- rep(prob, length(x))
    if(length(meanlog)==1) meanlog <- rep(meanlog, length(x))
    if(length(sdlog)==1) sdlog <- rep(sdlog, length(x))
    d <- 1-prob
    d[x>0] <- prob[x>0]*dlnorm(x[x>0], meanlog=meanlog[x>0], sdlog=sdlog[x>0])
    d
}
