orlm.forboot <- function(data, indices, index=index, ...){
    ## df.error chosen arbitrary, since it is irrelevant for b.restr
    dat <- data[indices,]
    bx <- orlm(cov.wt(dat[,-1],wt=dat$wt)$cov, df.error=10, index=index, ...)$b.restr
    ## intercept
    bi <- mean(dat[,2]) - sum(colMeans(dat[,3:ncol(dat)])*bx)
    return(c(bi,bx))
}
