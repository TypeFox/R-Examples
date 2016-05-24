qblnorm <-
function(p, prob, meanlog, sdlog)
{
    if(length(prob)==1) prob <- rep(prob, length(p))
    if(length(meanlog)==1) meanlog <- rep(meanlog, length(p))
    if(length(sdlog)==1) sdlog <- rep(sdlog, length(p))
    q <- rep(0, length(p))
    cases <- p > (1-prob)
    q[cases] <- qlnorm((prob[cases]+p[cases]-1)/prob[cases],
                       meanlog=meanlog[cases], sdlog=sdlog[cases])
    q
}
