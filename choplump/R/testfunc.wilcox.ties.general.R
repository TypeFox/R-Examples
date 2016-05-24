`testfunc.wilcox.ties.general` <-
function(d){
    r <- rank(d$W)
    n0 <- length(d$Z[d$Z==0])
    n1 <- length(d$Z[d$Z==1])
    STATISTIC <- sum(r[d$Z==0]) - n0 * (n0 + 1)/2
    NTIES <- table(r)
    SIGMA <- sqrt((n0 * n1/12) * ((n0 + n1 + 1) - 
            sum(NTIES^3 - NTIES)/((n0 + n1) * (n0 + n1 - 
              1))))
    out<-  (STATISTIC - n0 * n1/2)/SIGMA
    #pval<-pnorm(out)
    #return(pval)
    return(out) 
}

