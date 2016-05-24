"schoener.rtest" <- function(tr, keep, byburst=TRUE, nrep=500)
{
    ## Verifications
    if (!inherits(tr,"traj"))
        stop("tr should be of class traj")

    ## Removes the missing values
    tr <- tr[!is.na(tr$x),]
    tr <- tr[!is.na(tr$y),]

    ## Computes the observed schoener's ratio
    obs <- schoener(tr, keep, byburst)


    res <- list()
    trbis <- tr
    li <- split(tr[,c("x","y")], tr$id)
    if (byburst)
        li <- split(tr, tr$burst)

    foo <- function(trb) {
        x <- trb[,c("x","y")]
        r2 <- sum(((x[,1]-mean(x[,1]))^2) +
                  ((x[,2]-mean(x[,2]))^2))/(nrow(x) -1)
        d <- unclass(trb$date)
        diffd <- outer(d,d,"-")
        cons <- diffd>keep[1]&diffd<keep[2]
        m <- sum(cons)

        foobis <- function(i) {
            x <- trb[sample(1:nrow(trb), replace=FALSE),c("x","y")]
            t2tmp <- as.matrix(dist(x)^2)
            t2 <- sum(t2tmp[cons])/m
            rat <- t2/r2
            return(rat)
        }
        res <- unlist(lapply(1:nrep, foobis))
        return(res)
    }

    rr <- lapply(li, foo)
    rr <- lapply(1:length(rr),
                 function(x) as.randtest(rr[[x]], obs[x,1]))

    opar <- par(mfrow=n2mfrow(length(rr)))
    on.exit(par(opar))
    rr <- lapply(rr, function(x) { x$pvalue <- 1-x$pvalue; return(x)})
    names(rr) <- names(li)
    lapply(1:length(rr),
           function(x) plot(rr[[x]],
                            main = paste(names(rr)[x], ": P = ",
                            round(rr[[x]]$pvalue, 2), sep="")))
    return(rr)
}

