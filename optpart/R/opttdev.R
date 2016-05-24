opttdev <- function (veg,clustering,maxitr=100,minsiz=5) 
{
    clustering <- as.numeric(clustify(clustering))

    numplt <- nrow(veg)
    numspc <- ncol(veg)
    numcls <- length(table(clustering))
    sums <- rep(0,maxitr+1)
    numitr <- 0
    relsum <- rep(0,numcls)
    colsum <- rep(0,numcls)
    spcsum <- rep(0,numspc)
    tmpclu <- rep(0,numplt)
    tmp <- .Fortran('opttdev',
                    as.double(as.matrix(veg)),
                    as.integer(numplt),
                    as.integer(numspc),
                    clusid = as.integer(clustering),
                    as.integer(numcls),
                    as.integer(maxitr),
                    as.integer(minsiz),
                    sums = as.double(sums),
                    numitr = as.integer(numitr),
                    as.double(relsum),
                    as.double(colsum),
                    as.double(spcsum),
                    as.integer(tmpclu),
                    PACKAGE='optpart')
    out <- list()
    out$dev<- tmp$sums[1:tmp$numitr]
    out$numitr = tmp$numitr
    out$clustering <- tmp$clusid
    class(out) <- c('opttdev','clustering')
    return(out)
}
