exactpoissonPval<-function(x,T=1,r=1,relErr=1+1e-07,tsmethod="central"){
    if (length(r)>1 | length(x)>1 | length(T)>1) stop("x,T, and r should be of length 1")
    m<-r*T  
    if (tsmethod=="minlike"){
            ## copied from poisson.test
            if (m == 0) 
                PVAL<-(x == 0)
            else {
                d <- dpois(x, r * T)
                if (x == m) 
                   PVAL<-1
                else if (x < m) {
                  N <- ceiling(2 * m - x)
                  while (dpois(N, m) > d) N <- 2 * N
                  i <- seq.int(from = ceiling(m), to = N)
                  y <- sum(dpois(i, m) <= d * relErr)
                  PVAL<-ppois(x, m) + ppois(N - y, m, lower.tail = FALSE)
                }
                else {
                  i <- seq.int(from = 0, to = floor(m))
                  y <- sum(dpois(i, m) <= d * relErr)
                  PVAL<-ppois(y - 1, m) + ppois(x - 1, m, lower.tail = FALSE)
                }
            }
    } else if (tsmethod=="central"){
        PVAL<-min(1, 2*min(ppois(x, m), ppois(x - 
            1, m, lower.tail = FALSE)))
    } else if (tsmethod=="blaker"){
        lower.tail<-ppois(x,m)
        upper.tail<-ppois(x-1,m,lower.tail=FALSE)
        if (lower.tail>=upper.tail){
            ## find a new lower.tail that is as large as possible without being
            ## larger than upper.tail 
            LOWER.TAIL<-ppois(-1:x,m)
            xlower<-max((-1:x)[LOWER.TAIL<= upper.tail*relErr])
            PVAL<-upper.tail + ppois(xlower,m)
        } else {
            ## find a new upper.tail that is as large as possible without being
            ## larger than lower.tail
            ## Sept 21, 2012: fix for x=0 case, before multiplied by 2 and hung
            ## Now start at N=1
            if (x==0){ N<-1
            } else N<-2*x           
            while (ppois(N,m,lower.tail=FALSE)>lower.tail) N<-2*N
            y<-(x-1):N
            UPPER.TAIL<-ppois(y,m,lower.tail=FALSE)
            xupper<-min(y[UPPER.TAIL<=lower.tail*relErr]) 
            PVAL<-lower.tail+ ppois(xupper,m,lower.tail=FALSE)
        }
    }
    PVAL
}

exactpoissonPvals<-function(x,T=1,r=1,relErr=1+1e-07,tsmethod="central"){
    nr<-length(r)
    p<-rep(NA,nr)
    for (i in 1:nr){
        p[i]<-exactpoissonPval(x,T=T,r=r[i],relErr=relErr,tsmethod=tsmethod)
    }
    list(r=r,T=T,pvals=p)
}