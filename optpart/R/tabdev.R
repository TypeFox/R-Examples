tabdev <- function(x,...)
{
    UseMethod('tabdev')
}


tabdev.default <- function(x,clustering,nitr=999,...)
{
    if (!is.data.frame(x)) x <- data.frame(x)
    clustering <- as.integer(clustify(clustering))

    pltspc <- apply(x>0,2,sum)
    numplt <- nrow(x)
    numspc <- ncol(x)
    totdev <- 0.0
    spcdev <- rep(0.0,numspc)
    pval <- rep(0.0,numspc)
    ntypes <- max(clustering)
    relsum <- rep(0.0,ntypes)
    colsum <- rep(0.0,ntypes)
    spcsum <- rep(0.0,numspc)
    pclass <- rep(0,numplt)
    tclass <- rep(0,numplt)
    tmp <- .Fortran('tabdev',
        as.double(as.matrix(x)),
        as.integer(numplt),
        as.integer(numspc),
        as.integer(clustering),
        as.integer(ntypes),
        spcdev = as.double(spcdev),
        totdev = as.double(totdev),
        pval = as.double(pval),
        as.integer(nitr),
        as.double(relsum),
        as.double(colsum),
        as.double(spcsum),
        as.integer(pclass),
        as.integer(tclass),
        PACKAGE='optpart')
    tmp2 <- data.frame(names(x),pltspc,tmp$spcdev,round(tmp$pval,3))
    names(tmp2) <- c('species','numocc','deviance','pval')
    result <- list(spcdev=tmp2,totdev=tmp$totdev)
    class(result) <- 'tabdev'
    return(result)
}

tabdev.stride <- function(x,taxa,...)
{
    res <- rep(NA,ncol(x$clustering))
    for (i in 1:ncol(x$clustering)) {
        tmp <- tabdev(taxa,x$clustering[i])
        res[i] <- tmp$totdev
    }
    clusters <- x$seq
    tabdev <- res
    out <- data.frame(clusters,tabdev)
    out
}

summary.tabdev <- function (object,p=0.05,...) 
{
    if (class(object)!='tabdev') stop('You must pass an object of class objectdev')

    tmp <- object$spcdev
    tmp <- tmp[tmp$pval<=p,]
    tmp <- tmp[order(tmp$pval),]
    tmp
}

