### binning:  2014/06/10
# To manually bin the data

## either 'x' or 'freq' is given, not both/none
## if 'breaks' is not given, 'bw' should be given
## if 'from' is not given, take 'from = min(x)-runif(1)*bw'

## create a class "bdata":
###   plot(bdata) -- draw histogram
###   top.coded:  0=none, -1=left/lower,1=right/upper, 2=both
###   equalwidth: T/F
###   n: sample size
###   name: by 'x' or counts if 'x' is missing

## 2014/06/12: make binning a method to deal with 'histogram' object
## -- convert to 'bdata'.  Used in bde/lpsmooth

binning <- function(x, counts, nclass, breaks, bw)
    UseMethod("binning")

binning.histogram <-
    function(x, counts, nclass, breaks, bw)
{
    breaks <- x$breaks
    counts <- x$counts
    mids <- x$mids
    n <- sum(x$counts)
    nclass <- length(x$counts)
    tmp <- diff(x$breaks)
    if(any(tmp != tmp[1]))
        equalwidth <- FALSE
    else
        equalwidth <- TRUE

    if(!is.finite(breaks[1])){
        if(!is.finite(breaks[nclass+1])){
            top.coded <- 2
        }else{
            top.coded <- -1
        }
    }else{
        if(!is.finite(breaks[nclass+1])){
            top.coded <- 1
        }else{
            top.coded <- 0
        }
    }
    name <- x$xname
    
    structure(
        list(breaks = breaks,
             counts = counts,
             mids = mids,
             size = n,
             nclass = nclass,
             equalwidth=equalwidth,
             top.coded = top.coded,
             name = name,
             call = match.call()),
        class="bdata")
}

binning.default <- function(x, counts, nclass, breaks, bw)
{
    top.coded <- 0
    if(missing(x)){
        if(missing(counts))
            stop("Either 'x' or 'counts' should be given")
        name <- deparse(substitute(counts))
        counts <- counts[!is.na(counts)]
        stopifnot(is.numeric(counts))
        if(any(!is.finite(counts)))
            stop("Infinite count(s) not allowed")
        if(any(counts < 0))
            stop("Negative count(s) not allowed")
        nclass <- length(counts)
        if(counts[1]==0)
            warning("First class has frequency 0")
        if(counts[nclass]==0)
            warning("Last class has frequency 0")
        if(missing(breaks))
            stop("'breaks' (class limits) are missing")
        stopifnot(is.numeric(breaks))
        if(length(breaks)==1){
            stopifnot(is.finite(breaks))
            if(missing(bw))
                stop("Bin width 'bw' missing")
            stopifnot(is.numeric(bw))
            stopifnot(bw > 0)
            equalwidth <- TRUE
            breaks <- breaks + (0:nclass) * bw
            mids <- breaks[-1] - bw * 0.5
        }else{
            if(length(breaks) != nclass + 1)
                stop("'counts' and 'breaks' lengths not match")
            bws <- diff(breaks)
            if(any(bws <= 0))
                stop("Invalid 'breaks' value(s)")
            if(any(bws-bws[1] != 0)){
                equalwidth <- FALSE
            }else{
                equalwidth <- TRUE
            }
            mids <- breaks[-1] - bws * 0.5
            if(!is.finite(breaks[1])){
                if(!is.finite(breaks[nclass+1])){
                    top.coded <- 2
                }else{
                    top.coded <- -1
                }
            }else{
                if(!is.finite(breaks[nclass+1])){
                    top.coded <- 1
                }else{
                    top.coded <- 0
                }
            }
        }
        n <- sum(counts)
    }else{ #if 'x' not missing
        name <- deparse(substitute(x))
        equalwidth <- TRUE
        x <- x[!is.na(x)]
        if(!missing(counts)){
            counts <- counts[!is.na(counts)]
            stopifnot(is.numeric(counts))
            if(any(!is.finite(counts)))
                stop("Infinite count(s) not allowed")
            if(any(counts < 0))
                stop("Negative count(s) not allowed")
            if(any(counts - round(counts) != 0))
                stop("Non-integer 'counts' not allowed")
            if(length(x) != length(counts))
                stop("'x' and 'counts' have different lengths")
            x <- rep(x, counts)
        }

        x.finite <- is.finite(x)
        if(any(!x.finite))
            warning(paste(sum(!x.finite), " infinite value(s) found and removed"))
        x <- x[x.finite]
        ## sample size
        n <- length(x)
        
        if(missing(breaks)){
            if(missing(bw)){
                if(missing(nclass))
                    stop("One of 'nclass', 'breaks' and 'bw' needs to be specified")
                stopifnot(is.numeric(nclass))
                stopifnot(length(nclass) == 1)
                stopifnot(nclass > 0)
                nclass <- round(nclass)
                bw2 <- diff(range(x))/nclass
                x0 <- min(x) - runif(1)*bw2
                x1 <- max(x) + runif(1)*bw2
                bw <- (x1-x0)/nclass
                breaks <- seq(x0, x1, length=nclass+1)
                mids <- breaks[-1] - bw * 0.5
            }else{
                stopifnot(is.numeric(bw))
                stopifnot(length(bw) == 1)
                stopifnot(bw>0)
                if(!missing(nclass))
                    warning("'nclass' not used")
                x0 <- min(x) - runif(1)*bw
                nclass <- ceiling((max(x)-x0)/bw)
                equalwidth <- TRUE
                breaks <- x0 + (0:nclass)*bw
                mids <- breaks[-1] - bw * 0.5
            }
        }else{ # breaks specified
            stopifnot(is.numeric(breaks))
            if(length(breaks)==1){
                stopifnot(is.finite(breaks))
                x0 <- breaks
                if(x0 > min(x))
                    stop("Small 'x' value(s) not in the 1st class")
                if(missing(bw)){
                    if(missing(nclass))
                        stop("None of 'nclass', 'breaks' and 'bw' specified")
                    stopifnot(is.numeric(nclass))
                    stopifnot(length(nclass) == 1)
                    stopifnot(nclass > 0)
                    nclass <- round(nclass)
                    bw2 <- diff(range(x))/nclass
                    x1 <- max(x) + runif(1)*bw2
                    bw <- (x1-x0)/nclass
                    equalwidth <- TRUE
                    breaks <- seq(x0, x1, length=nclass+1)
                    mids <- breaks[-1] - bw * 0.5
                }else{
                    stopifnot(is.numeric(bw))
                    stopifnot(length(bw) == 1)
                    stopifnot(bw>0)
                    if(!missing(nclass))
                        warning("'nclass' not used")
                    nclass <- ceiling((max(x)-x0)/bw)
                    equalwidth <- TRUE
                    breaks <- x0 + (0:nclass)*bw
                    mids <- breaks[-1] - bw * 0.5
                }
            }else{
                nclass <- length(breaks)-1
                bws <- diff(breaks)
                if(any(bws <= 0))
                    stop("Invalid 'breaks' value(s)")
                if(any(bws-bws[1] != 0)){
                    equalwidth <- FALSE
                }else{
                    equalwidth <- TRUE
                }
                mids <- breaks[-1] - bws * 0.5
                if(!is.finite(breaks[1])){
                    if(!is.finite(breaks[nclass+1])){
                        top.coded <- 2
                    }else{
                        top.coded <- -1
                    }
                }else{
                    if(!is.finite(breaks[nclass+1])){
                        top.coded <- 1
                    }else{
                        top.coded <- 0
                    }
                }
            }
        }
        res <- hist(x, breaks=breaks,plot=FALSE)
        counts <- res$counts
    }
    structure(
        list(breaks = breaks,
             counts = counts,
             mids = mids,
             size = n,
             nclass = nclass,
             equalwidth=equalwidth,
             top.coded = top.coded,
             name = name,
             call = match.call()),
        class="bdata")  
}

.summary.bdata <- function(x){
    if(class(x)=="histogram")
        x <- binning(x)
    stopifnot(class(x)=="bdata")
    ## compute exact quantiles
    Fn <- cumsum(x$counts)/sum(x$counts)
    qexact <- .qtlexact(Fn=Fn, breaks=x$breaks)
    
    ## compute approximate quantiles
    qlevels <- c(0.1,0.25,0.5,0.75,0.9)
    nclass <- x$nclass
    if(!is.finite(x$breaks[1])){
        q.small <- qlevels < Fn[1]
        if(any(q.small)){
            j <- sum(q.small)
            tmp <- seq(Fn[1], Fn[2], length=j+1)
            qlevels[1:j] <- tmp[-(j+1)]
        }
    }

    if(!is.finite(x$breaks[nclass+1])){
        q.small <- qlevels > Fn[nclass-1]
        if(any(q.small)){
            j <- sum(q.small)
            tmp <- seq(Fn[nclass-2], Fn[nclass-1], length=j+1)
            qlevels[(6-j):5] <- tmp[-1]
        }
    }

    qtls1 <- lapply(qlevels,.qtlapprox,Fn=Fn,breaks=x$breaks)
    qappr <- list(levels=qlevels,qtls=as.numeric(qtls1))


    ## compute mean and variance
    ## moments
    breaks <- x$breaks
    counts <- x$counts
    k <- x$nclass
    if(x$top.coded==2){
        breaks <- breaks[-c(1,k+1)]
        counts <- counts[-c(1,k)]
    }else if(x$top.coded==1){
        breaks <- breaks[-(k+1)]
        counts <- counts[-(k)]
    }else if(x$top.coded==-1){
        breaks <- breaks[-1]
        counts <- counts[-1]
    }
    k <- length(breaks)
    mids <- 0.5*(breaks[-1]+breaks[-k])
    mu <- sum(mids * counts)/sum(counts)
    M2 <- sum((mids-mu)^2*counts)/(sum(counts)-1)
    ## qtls1 are approximated, qtls2 are exact (if exist)
    res <- list(mean=mu, var=M2, median=qappr$qtls[3],
                exact=qexact, approx=qappr)
}

print.bdata <- function(x,...){
    cat("Call:  ", deparse(x$call), sep = "")
    cat("\n  Data: ", x$name, "\t(n=", x$size, 
        "; nclass=", x$nclass, ")\n", sep="")
    cat("  Finite limits: ", x$top.coded==0,
        "; equal bw: ", x$equalwidth,
        "\n", sep="")
    res <- .summary.bdata(x)
    cat("  Mean = ", res$mean, "\tMedian = ", res$median,
        "\tSD = ", round(sqrt(res$var),3),
        ")\n", sep = "")
    cat("  Quantiles (approximate):\n")
    tmp <- res$approx$qtls
    names(tmp) <- round(res$approx$levels,4)
    print(tmp)
    
    cat("  Quantiles (exact):\n")
    tmp <- res$exact$qtls
    names(tmp) <- round(res$exact$levels, 4)
    print(tmp)

    lbs <- paste("(", x$breaks[-(x$nclass+1)],",  ",
                 x$breaks[-1], "]", sep="")
    n <- length(lbs)
    if(x$top.coded!=-1){
        l <- nchar(lbs[1])
        lbs[1] <- paste("[", substr(lbs[1],2,l),sep='')
    }
    if(x$top.coded==1||x$top.coded==2){
        l <- nchar(lbs[n])
        lbs[n] <- paste(substr(lbs[n],1,l-1),")",sep='')
    }
    perc <- x$counts/sum(x$counts)*100
    cperc <- cumsum(x$counts)/sum(x$counts)*100
    tmp <- data.frame("Classes"=lbs,
                      "Frequencies"=x$counts,
                      "Perc"=round(perc,2),
                      "CumPerc"=round(cperc,2))
    print(tmp)
    cat("\n")
}

.bdataTohist <- function(x){
    if(x$top.coded != 0)
        x$name <- paste(x$name, " (top-coded)", sep='')
    counts <- x$counts
    breaks <- x$breaks
    bw <- diff(breaks)
    nclass <- length(counts)
    
    if(x$top.coded == -1 || x$top.coded == 2){
        bw1 <- 5 * counts[1] * bw[2] / counts[2]
        breaks[1] <- breaks[2] - bw1
        bw <- diff(breaks)
    }
        
    if(x$top.coded == 1 || x$top.coded == 2){
        bw1 <- 5 * counts[nclass]*bw[nclass-1]/counts[nclass-1]
        breaks[nclass+1] <- breaks[nclass] + bw1
        bw <- diff(breaks)
    }
    mids <- 0.5*(breaks[-1]+breaks[-(nclass+1)])

    res <- structure(
        list(breaks = breaks,
             counts = counts,
             density = counts/bw/sum(counts),
             mids = mids,
             xname = x$name,
             equidist=FALSE),
        class="histogram")
}

plot.bdata <- function(x,...){
    res <- .bdataTohist(x)
    plot(res,...)
}

 #> ofc = rnorm(1000,34.5,1.5)
#> GLD.fit.default(ofc)
#           [,1]    [,2]     [,3]     [,4]     [,5]
#qlevels  0.1000  0.2500  0.50000  0.75000  0.90000
#qtls    32.5443 33.5065 34.52928 35.51963 36.42914

.qtlapprox <- function(q, Fn, breaks){
    K <- length(breaks)
    if(q<0 || q>1){
        res <- NA
        warning("Invalid quantile level")
    }else if(q==0){
        res <- breaks[1]
    }else if(q==1){
        res <- breaks[K]
    }else{
        if(any(Fn == q)){
            i <- which(Fn == q)[1]
            res <- breaks[i+1]
        }else{
            i <- which(Fn > q)[1]
            if(i <= 1){
                if(!is.finite(breaks[1]))
                    res <- NA
                else
                    res <- breaks[1] + (breaks[2]-breaks[1])*q/Fn[1]
            }else if(i > K - 2){
                if(!is.finite(breaks[K]))
                    res <- NA
                else
                    res <- breaks[K-1]+(breaks[K]-breaks[K-1])*(q-Fn[K-2])/(1-Fn[K-2])
            }else
                res <- breaks[i] + (breaks[i+1]-breaks[i])*(q-Fn[i-1])/(Fn[i]-Fn[i-1])
        }
    }
    res
}

.qtlexact <- function(Fn, breaks){
    ##    qlevels <- c(0.1,0.25,0.5,0.75,0.9)
    Fn <- c(0, Fn)
    K <- length(breaks)
    if(K <= 6){
        res <- NA; qlevels <- NA
        x.finite <- is.finite(breaks)
        res <- list(levels=Fn[x.finite], qtls=breaks[x.finite])
    }else{
        ## find the median q3 among the breaks points from 4=(1+3)...(K-3)
        qdiff <- abs(Fn - 0.5)
        i3 <- which(qdiff == min(qdiff))[1]
        if(i3 < 4) i3 <- 4
        else if(i3 > K-3) i3 <- K-3
        q3 <- Fn[i3]; x3 <- breaks[i3]
        ## find the median q2 among the breaks points from 3=(1+2)...(i3-1)
        qdiff <- abs(Fn - 0.25)
        i2 <- which(qdiff == min(qdiff))[1]
        if(i2 < 3) i2 <- 3
        else if(i2 > i3-1) i2 <- i3-1
        q2 <- Fn[i2]; x2 <- breaks[i2]
        ## find the median q1 among the breaks points from 2=(1+1)...(i2-1)
        qdiff <- abs(Fn - 0.10)
        i1 <- which(qdiff == min(qdiff))[1]
        if(i1 < 2) i1 <- 2
        else if(i1 > i2-1) i1 <- i2-1
        q1 <- Fn[i1]; x1 <- breaks[i1]
        ## find the median q4 among the breaks points from (i3+1)...(K-2)
        qdiff <- abs(Fn - 0.75)
        i4 <- which(qdiff == min(qdiff))[1]
        if(i4 < i3+1) i4 <- i3+1
        else if(i4 > K-2) i4 <- K-2
        q4 <- Fn[i4]; x4 <- breaks[i4]
        ## find the median q5 among the breaks points from (i4+1)...(K-1)
        qdiff <- abs(Fn - 0.90)
        i5 <- which(qdiff == min(qdiff))[1]
        if(i5 < i4+1) i5 <- i4+1
        else if(i5 > K-1) i5 <- K-1
        q5 <- Fn[i5]; x5 <- breaks[i5]
        res <- list(levels=c(q1,q2,q3,q4,q5), qtls=c(x1,x2,x3,x4,x5))
    }
    res
}

.compact <- function(x){
    x0 <- x
    x <- x0$breaks; y <- x0$counts
    while(any(y <  5)){
        ymin <- min(y); ny <- length(y)
        imin <- which(y == ymin)[1]
        if(imin == 1){
            x <- x[-1]; y <- y[-1]
            y[1] <- y[1] + ymin
        }else if(imin == ny){
            x <- x[-(ny+1)]; y <- y[-ny]
            y[ny-1] <- y[ny-1] + ymin
        }else{
            if(y[imin-1] < y[imin+1]){
                x <- x[-imin]; y <- y[-imin]
                y[imin-1] <- y[imin-1] + ymin
            }else{
                x <- x[-(imin+1)]; y <- y[-imin]
                y[imin] <- y[imin] + ymin
            }
        }
    }
    x0$counts <- y; x0$breaks <- x
    nx <- length(x)
    x0$nclass <- nx-1
    x0$mids <- (x[-1]+x[-nx])*0.5
    x0
}
