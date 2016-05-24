## "myplotDD" is a revised version of ".myddplot" in "plot-utils.R" in the
## package "rrcov". In "myplotDD", id.n and ind are printed out.
## 
## id.n : Number of observations to identify by a label. 
##        If not supplied, the number of observations with robust distance 
##        larger than cutoff is used.
## ind  : The index of rd whose values are larger than cutoff.
## 
## 
myplotDD <- function(x, cutoff, id.n) {
    ##  Distance-Distance Plot:
    ##  Plot the vector y=rd (robust distances) against
    ##  x=md (mahalanobis distances). Identify by a label the id.n
    ##  observations with largest rd. If id.n is not supplied, calculate
    ##  it as the number of observations larger than cutoff. Use cutoff
    ##  to draw a horisontal and a vertical line. Draw also a dotted line
    ##  with a slope 1.

    data <- getData(x)
    ##  parameters and preconditions
    if(is.vector(data) || is.matrix(data)) {
        if(!is.numeric(data))
            stop( "x is not a numeric dataframe or matrix.")
    } else if(is.data.frame(data)) {
        if(!all(sapply(data,data.class) == "numeric"))
            stop( "x is not a numeric dataframe or matrix.")
    }

    n <- dim(data)[1]
    p <- dim(data)[2]

    if(length(getCenter(x))  == 0 ||  length(getCov(x)) == 0)
        stop( "Invalid object: attributes center and cov missing!")

    if(length(getCenter(x))  != p)
        stop( "Data set and provided center have different dimensions!")

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))
        cat("cutoff =\n"); show(cutoff)

    if(!missing(id.n) && !is.null(id.n)) {
        id.n <- as.integer(id.n)
        if(id.n < 0 || id.n > n)
            stop(sQuote("id.n")," must be in {1,..,",n,"}")
    }

    ccov <- CovClassic(data)
    md <- rd <- NULL
    md <- sqrt(getDistance(ccov))
    rd <- sqrt(getDistance(x))

# .myddplot(md, rd, cutoff=cutoff, id.n=id.n) # distance-distance plot
    n <- length(md)
    if(missing(id.n)){
        id.n <- length(which(rd>cutoff)); 
        cat("id.n <- length(which(rd>cutoff))\n")
        cat("id.n =", id.n, "\n")
    }

    xlab <- "Mahalanobis distance"
    ylab <- "Robust distance"
    plot(md, rd, xlab=xlab, ylab=ylab, type="p")

#    .label(md,rd,id.n)
    x=md; y=rd; 
    cat("Here y is the robust distance (rd).\n");
    xrange <- par("usr")
    xrange <- xrange[2] - xrange[1]
    n <- length(y)
    sort.y <- sort(y, index.return=TRUE); cat("sort.y =\n"); show(sort.y)
    ind.sort.y <- sort.y$ix; 
    ind <- ind.sort.y[(n-id.n+1):n]; cat("ind =\n"); show(ind)
    text(x[ind] + xrange/50, y[ind], ind)

    abline(0, 1, lty=2)
    abline(v=cutoff)
    abline(h=cutoff)
    title(main="My Distance-Distance Plot")

    result=list(cutoff = cutoff, id.n = id.n, sort.y = sort.y, ind = ind)
    result
}

