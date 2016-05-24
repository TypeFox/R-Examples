envrtest <- function (set,env,numitr=1000,minval=0, replace=FALSE,
                plotit=TRUE, main=paste(deparse(substitute(set)),' on ',
                deparse(substitute(env)))) 
{
    prob <- 0
    if (!is.logical(set)) {
        mask <- (set > minval)
    } else {
        mask <- set
    }
    omin <- min(env[mask])
    omax <- max(env[mask])
    odiff <- omax - omin
    sdiff <- rep(0,numitr-1)
    for (i in 1:numitr-1) {
        tmp <- sample(1:length(env),sum(mask),replace=replace)
        nullmin <- min(env[tmp])
        nullmax <- max(env[tmp])
        null <- nullmax - nullmin
        if (null <= odiff) prob <- prob + 1
        sdiff[i] <- null
    }
    prob <- min(1,(prob+1)/numitr)
    if (plotit) {
        plot(sort(sdiff),ylim=c(min(odiff,sdiff),max(sdiff)),
             ylab="Within-Set Difference",
             main=main)
        abline(odiff,0,col=2)
        text(0,max(sdiff),paste("p = ",format(prob,digits=3)),adj=0)
    }
    out <- list()
    out$diff <- odiff
    out$prob <- prob
    attr(out,'call') <- match.call()
    invisible(out)
}
