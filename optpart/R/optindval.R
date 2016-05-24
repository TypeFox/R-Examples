optindval <- function (veg,clustering,maxitr=100,minsiz=5) 
{
    numplt <- nrow(veg)
    numspc <- ncol(veg)
    clustering <- as.integer(clustify(clustering))

    numcls <- max(clustering)
    relfrq <- matrix(0,nrow=numspc,ncol=numcls)
    relabu <- matrix(0,nrow=numspc,ncol=numcls)
    indval <- matrix(0,nrow=numspc,ncol=numcls)
    pval <- rep(0,numspc)
    indcls <- rep(0,numspc)
    clstab <- rep(0,numcls)
    maxcls <- rep(0,numspc)
    sumind <- 0
    sums <- rep(0,maxitr+1)
    numitr <- 0
    tmpclu <- rep(0,numplt)
    res <- .Fortran('optindval',
                   as.double(as.matrix(veg)),
                   as.integer(numplt),
                   as.integer(numspc),
                   clustering = as.integer(clustering),
                   as.integer(numcls),
                   as.double(relfrq),
                   as.double(relabu),
                   as.double(indval),
                   as.double(indcls),
                   as.integer(clstab),
                   as.integer(maxcls),
                   as.double(sumind),
                   as.integer(maxitr),
                   as.integer(minsiz),
                   sums = as.double(sums),
                   numitr = as.integer(numitr),
                   as.integer(tmpclu),
                   PACKAGE='optpart')
    out <- indval(veg,res$clustering)
    out$clustering <- res$clustering
    out$sums <- res$sums[1:res$numitr]
    out$numitr <- res$numitr
    class(out) <- c('optindval','clustering')
    return(out)
}
