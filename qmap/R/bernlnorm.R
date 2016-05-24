dbernlnorm <- function(x, prob, meanlog, sdlog){
    if(length(prob)==1) prob <- rep(prob, length(x))
    if(length(meanlog)==1) meanlog <- rep(meanlog, length(x))
    if(length(sdlog)==1) sdlog <- rep(sdlog, length(x))
    d <- 1-prob
    d[x>0] <- prob[x>0]*dlnorm(x[x>0], meanlog=meanlog[x>0], sdlog=sdlog[x>0])
    d
}


pbernlnorm <- function(q, prob, meanlog, sdlog){
    if(length(prob)==1) prob <- rep(prob, length(q))
    if(length(meanlog)==1) meanlog <- rep(meanlog, length(q))
    if(length(sdlog)==1) sdlog <- rep(sdlog, length(q))
    p <- 1-prob
    p[q>0] <- 1-prob[q>0]+prob[q>0]*plnorm(q[q>0], meanlog=meanlog[q>0],
                                           sdlog=sdlog[q>0])
    p 
}

qbernlnorm <- function(p, prob, meanlog, sdlog){
    if(length(prob)==1) prob <- rep(prob, length(p))
    if(length(meanlog)==1) meanlog <- rep(meanlog, length(p))
    if(length(sdlog)==1) sdlog <- rep(sdlog, length(p))
    q <- rep(0, length(p))
    cases <- p > (1-prob)
    q[cases] <- qlnorm((prob[cases]+p[cases]-1)/prob[cases],
                       meanlog=meanlog[cases], sdlog=sdlog[cases])
    q
}

rbernlnorm <- function(n, prob, meanlog, sdlog){
    if(max(length(prob), length(meanlog), length(sdlog)) > 1)
        stop("parameters must be of length 1")
    p <- runif(n)
    q <- rep(0, length(p))
    cases <- p > (1-prob)
    q[cases] <- rlnorm(sum(cases), meanlog=meanlog, sdlog=sdlog)
    q
}

