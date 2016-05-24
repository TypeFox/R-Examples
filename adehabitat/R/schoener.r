"schoener" <- function(tr, keep, byburst=TRUE)
{
    ## Verifications
    if (!inherits(tr, "traj"))
        stop("tr should be of class traj")

    ## Remove the missing values
    tr <- tr[!is.na(tr$x),]
    tr <- tr[!is.na(tr$y),]

    ## splits per burst or id
    li <- split(tr, tr$id)
    if (byburst)
        li <- split(tr, tr$burst)

    ## This function computes the schoener ratio
    ## for each element of this list
    foo <- function(tr) {

        ## Computation of r2
        d <- unclass(tr$date)
        x <- tr[,c("x","y")]
        r2 <- sum(((x[,1]-mean(x[,1]))^2) +
                  ((x[,2]-mean(x[,2]))^2))/(nrow(x) -1)

        ## Computation of t2
        diffd <- outer(d,d,"-")
        t2tmp <- as.matrix(dist(x)^2)
        cons <- diffd>keep[1]&diffd<keep[2]
        t2 <- sum(t2tmp[cons])/sum(cons)

        ## The ratio
        rat <- t2/r2
        n <- nrow(x)
        m <- sum(cons)

        ## Output
        return(c(rat, n, m))
    }

    ## The function is applied to each element of the list
    rr <- do.call("rbind", lapply(li, foo))
    rr <- as.data.frame(rr)
    row.names(rr) <- names(li)
    names(rr) <- c("value","n","m")

    ## Output
    return(rr)
}

