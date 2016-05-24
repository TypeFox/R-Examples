dbernexp <- function(x, prob, rate){
    if(length(prob)==1) prob <- rep(prob, length(x))
    if(length(rate)==1) rate <- rep(rate, length(x))   
    d <- 1-prob
    d[x>0] <- prob[x>0]*dexp(x[x>0], rate=rate[x>0])
    d
}


pbernexp <- function(q, prob, rate){
    if(length(prob)==1) prob <- rep(prob, length(q))
    if(length(rate)==1) rate <- rep(rate, length(q))
    p <- 1-prob
    p[q>0] <- 1-prob[q>0]+prob[q>0]*pexp(q[q>0],
                                         rate=rate[q>0])
    p 
}

qbernexp <- function(p, prob, rate){
    if(length(prob)==1) prob <- rep(prob, length(p))
    if(length(rate)==1) rate <- rep(rate, length(p))
    q <- rep(0, length(p))
    cases <- p > (1-prob)
    q[cases] <- qexp((prob[cases]+p[cases]-1)/prob[cases],
                     rate=rate[cases])
    q
  }

rbernexp <- function(n, prob, rate){
    if(max(length(prob), length(rate)) > 1)
        stop("parameters must be of length 1")
    p <- runif(n)
    q <- rep(0, length(p))
    cases <- p > (1-prob)
    q[cases] <- rexp(sum(cases), rate=rate)
    q
}


