optsil <- function(x, dist, maxitr)
{
    UseMethod('optsil')
}

opts.core <- function (dist, clustering, maxitr = 100) 
{
    if (class(dist) != 'dist') 
        stop('You must pass an object of classs dist as the first argument')
    if (max(dist) > 1) dist <- dist/max(dist)
    sim <- 1 - as.matrix(dist)
    clustering <- as.numeric(clustify(clustering))
    numplt <- length(clustering)
    numclu <- max(clustering)
    sils <- rep(0, maxitr)
    numitr <- 0
    pltsil <- rep(0,numplt)
    tmpclu <- rep(0,numplt)
    simptc <- matrix(0,nrow=numplt,ncol=numclu) 
    nabor <- rep(0,numplt)
    sumnum <- rep(0,numclu)
    sumden <- rep(0,numclu)
    res <- .Fortran('optsil', 
               as.double(sim), 
               clustering = as.integer(clustering), 
               as.integer(numplt), 
               as.integer(numclu), 
               as.integer(maxitr), 
               sils = as.double(sils), 
               numitr = as.integer(numitr), 
               as.double(simptc),
               as.double(pltsil),
               as.integer(tmpclu),
               as.integer(nabor),
               as.double(sumnum),
               as.integer(sumden),
               PACKAGE = 'optpart')
    out <- list()
    out$clustering <- as.numeric(factor(res$clustering))
    out$sils <- (res$sils/numplt)[1:res$numitr]
    out$numitr <- res$numitr
    class(out) <- c('optsil', 'clustering')
    out
}

optsil.default <- function(x,dist,maxitr=100)
{
    if (class(dist) != 'dist') {
        stop('You must pass an object of class dist')
    }
    
    if (is.numeric(x) && length(x) == 1) {
        out <- opts.core(dist,sample(1:x,attr(dist,'Size'),
            replace=TRUE),maxitr)
    } else {
        clustering <- as.numeric(clustify(x))
        out <- opts.core(dist,clustering,maxitr)
    }

    attr(out,'class') <- c('optsil','clustering')
    attr(out,'call') <- match.call()
    out
}

optsil.stride <- function(x,dist,maxitr=100)
{
    if (class(x) != 'stride') 
        stop('You must pass an object of class stride')
    res <- matrix(NA,nrow=nrow(x$clustering),ncol=ncol(x$clustering))
    for (i in 1:ncol(x$clustering)) {
        tmp <- opts.core(dist,x$clustering[,i],maxitr)
        res[,i] <- tmp$clustering
    }
    out <- data.frame(res)
    names(out) <- as.character(x$seq)
    row.names(out) <- row.names(x$clustering)
    out <- list(clustering=out,seq=x$seq)
    attr(out,'class') <- 'stride'
    attr(out,'call') <- match.call()
    out
}
