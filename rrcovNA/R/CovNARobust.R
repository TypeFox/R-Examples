##  control can be a character specifying the name of the estimate, one of:
##  auto, mcd, ogk, sde, sfast, surreal, bisquare, rocke
##  If no control object is given or 'auto' is selected, the choice of the
##  estimator will depend on the size of the data: see CovRobust() in package rrcov.
##
CovNARobust <- function(x, control, impMeth=c("norm" , "seq", "rseq"))
{
    impMeth <- match.arg(impMeth)

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if(!is.matrix(x))
        x <- matrix(x, length(x), 1,
            dimnames = list(names(x), deparse(substitute(x))))

    xcall <- match.call()
    n <- nrow(x)
    p <- ncol(x)

    ## drop all rows which contain only missings
    na.x <- rowSums(ifelse(is.na(x),1,0)) == ncol(x)
    ok <- !na.x
    x <- x[ok, , drop = FALSE]

    dimn <- dimnames(x)
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    ximp <- .imputation(x, impMeth = impMeth)
    CovRobust(ximp, control)
}

setMethod("summary", "CovNARobust", function(object, ...){

    new("SummaryCovNARobust", covobj=object, evals=eigen(object@cov)$values)

})


setMethod("show", "SummaryCovNARobust", function(object){

    cat("\nCall:\n")
    print(object@covobj@call)

    digits = max(3, getOption("digits") - 3)
    cat("\nRobust Estimate of Location: \n")
    print.default(format(getCenter(object), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nRobust Estimate of Covariance: \n")
    print.default(format(getCov(object), digits = digits), print.gap = 2, quote = FALSE)

    cat("\nEigenvalues of covariance matrix: \n")
    print.default(format(getEvals(object), digits = digits), print.gap = 2, quote = FALSE)

    cat("\nRobust Distances: \n")
    dd <- getDistance(object)
    print.default(format(as.vector(dd$d), digits = digits), print.gap = 2, quote = FALSE)
})

setMethod("plot", signature(x="CovNARobust", y="missing"), function(x, y="missing",
                                which=c("dd", "distance", "qqchi2", "tolEllipsePlot", "screeplot", "pairs", "all", "xydistance", "xyqqchi2"),
                                classic= FALSE,
                                ask = (which=="all" && dev.interactive(TRUE)),
                                cutoff,
                                id.n,
                                tol = 1e-7, ...)
{
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

    ## Check for singularity of the cov matrix
    if(isSingular(x))
        stop("The covariance matrix is singular!")

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(!missing(id.n) && !is.null(id.n)) {
        id.n <- as.integer(id.n)
        if(id.n < 0 || id.n > n)
            stop(sQuote("id.n")," must be in {1,..,",n,"}")
    }

    ccov <- CovNAClassic(data)
##    print(getDistance(ccov))

    md <- rd <- NULL
    if(!isSingular(ccov))
        md <- sqrt(getDistance(ccov)$d)
    if(!isSingular(x))
        rd <- sqrt(getDistance(x)$d)

    which <- match.arg(which)
    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))

    mymiss <- .xmiss(x@X)

    ## distance-distance plot: here we need both robust and mahalanobis distances
    if((which == "all" || which == "dd") && !is.null(md) && !is.null(rd)) {
        .myddplotNA(md, rd, mymiss, cutoff=cutoff, id.n=id.n) # distance-distance plot
    }

    ## index plot of mahalanobis distances
    if((which == "all" || which == "distance")  && !is.null(rd)) {
        ylim <- NULL
        if(classic && !is.null(md)) {
            opr <- if(prod(par("mfrow")) == 1) par(mfrow=c(1,2), pty="m") else list()

            ##VT::10.11.2007 - set same scale on both plots
            ylim <- c(min(rd,md), max(md,rd))
        }

        .mydistplotNA(rd, mymiss, cutoff, id.n=id.n, ylim=ylim)                   # index plot of robust distances
        if(classic && !is.null(md)) {
            .mydistplotNA(md, mymiss, cutoff, classic=TRUE, id.n=id.n, ylim=ylim) # index plot of mahalanobis distances
            par(opr)
        }
    }

    ## lattice: index plot of mahalanobis distances
    if(which == "xydistance"  && !is.null(rd)) {
        print(.xydistplotNA(x, cutoff=cutoff, ...))      # lattice:  index plot of robust distances

    }

    ## qq-plot of the mahalanobis distances versus the
    ## quantiles of the chi-squared distribution
    if((which == "all" || which == "qqchi2")  && !is.null(rd)) {
        if(classic && !is.null(md)) {
            opr <- if(prod(par("mfrow")) == 1) par(mfrow=c(1,2), pty="m") else list()
        }
        .qqplotNA(rd, p, mymiss, cutoff=cutoff, id.n=id.n) # qq-plot of the robust distances versus the
                                                 # quantiles of the chi-squared distribution
        if(classic && !is.null(md)) {
            .qqplotNA(md, p, mymiss, cutoff=cutoff, classic=TRUE, id.n=id.n)
                                                 # qq-plot of the mahalanobis distances
            par(opr)
        }
    }

    ## lattice: qq-plot of the mahalanobis distances versus the
    ##          quantiles of the chi-squared distribution
    if(which == "xyqqchi2"  && !is.null(rd)) {
        print(.xyqqchi2NA(x, cutoff=cutoff, ...))        # lattice:  qq-plot of the distances versus
    }


    if(which == "tolEllipsePlot" || which == "pairs") {
        if(which == "tolEllipsePlot" & length(dim(data)) >= 2 && dim(data)[2] == 2){
            if(!is.null(rd)){
                if(classic &&  !is.null(md))
                    .tolellipseNA(rcov=x, ccov = ccov, cutoff=cutoff, id.n=id.n, tol=tol, ...)
                else
                    .tolellipseNA(rcov=x, cutoff=cutoff, id.n=id.n, tol=tol, ...)
            }
        }else if(length(dim(data)) >= 2 && dim(data)[2] <= 10)
        {
            .rrpairsNA(x, ...)
        }else if(which != "all")
            warning("Warning: For tolerance ellipses the dimension must be less than 10!")
    }

    if(which == "all" || which == "screeplot") {
        myscreeplot(ccov=ccov, rcov=x)
    }
}) ## end { plot("CovRobust") }

.labelNA <- function(x, y, xmiss, id.n=3) {
    if(id.n > 0) {
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        n <- length(y)
        ind <- sort(y, index.return=TRUE)$ix
        ind <- ind[(n-id.n+1):n]
        text(x[ind] + xrange/50, y[ind], ind, col=ifelse(xmiss[ind],"red","black"))
    }
}

.myddplotNA <- function(md, rd, xmiss, cutoff, id.n) {
    ##  Distance-Distance Plot:
    ##  Plot the vector y=rd (robust distances) against
    ##  x=md (mahalanobis distances). Identify by a label the id.n
    ##  observations with largest rd. If id.n is not supplied, calculate
    ##  it as the number of observations larger than cutoff. Use cutoff
    ##  to draw a horisontal and a vertical line. Draw also a dotted line
    ##  with a slope 1.
    n <- length(md)
    if(missing(id.n))
        id.n <- length(which(rd>cutoff))
    xlab <- "Mahalanobis distance"
    ylab <- "Robust distance"
    plot(md, rd, xlab=xlab, ylab=ylab, type="n")
    points(md, rd, type="p", col=ifelse(xmiss, "red", "black"))
   .labelNA(md, rd, xmiss, id.n)
    abline(0, 1, lty=2)
    abline(v=cutoff)
    abline(h=cutoff)

    title(main="Distance-Distance Plot")
}

.mydistplotNA <- function(x, xmiss, cutoff, classic = FALSE, id.n, ylim=NULL) {
    ## VT::10.11.2007 - change "Squared Robust distance" to "Robust distance"
    ## VT::10.11.2007 - Add parameter ylim to make possible robust and
    ##  classical plot to use the same y-scale

    ##  Index Plot:
    ##  Plot the vector x (robust or mahalanobis distances) against
    ##  the observation indexes. Identify by a label the id.n
    ##  observations with largest value of x. If id.n is not supplied,
    ##  calculate it as the number of observations larger than cutoff.
    ##  Use cutoff to draw a horisontal line.
    ##  Use classic=FALSE/TRUE to choose the label of the vertical axes

    n <- length(x)
    if(missing(id.n))
        id.n <- length(which(x>cutoff))

    ylab <- paste(if(classic) "Mahalanobis" else "Robust", "distance")
    plot(x, ylab=ylab, xlab="Index", type="n", ylim=ylim)
    points(x, type="p", col=ifelse(xmiss,"red","black"))
    .labelNA(1:n, x, xmiss, id.n)
    abline(h=cutoff)

    title(main="Distance Plot")
}

.qqplotNA <- function(x, p, xmiss, cutoff, classic=FALSE, id.n) {
    ##  Chisquare QQ-Plot:
    ##  Plot the vector x (robust or mahalanobis distances) against
    ##  the square root of the quantiles of the chi-squared distribution
    ##  with p degrees of freedom.
    ##  Identify by a label the id.n observations with largest value of x.
    ##  If id.n is not supplied, calculate it as the number of observations
    ##  larger than cutoff.
    ##  Use classic=FALSE/TRUE to choose the label of the vertical axes


    ##  parameters and preconditions

    n <- length(x)

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(missing(id.n))
        id.n <- length(which(x>cutoff))

    qq <- sqrt(qchisq(((1:n)-1/3)/(n+1/3), p))

    x <- sort(x, index.return=TRUE)
    ix <- x$ix
    x <- x$x

    if(classic)
        ylab="Mahalanobis distance"
    else
        ylab="Robust distance"

    ## xlab <- "Square root of the quantiles of the chi-squared distribution"
    xlab <- eval(substitute(expression(paste("Sqrt of the quantiles of the ", chi^2, " distribution"))))
    plot(qq, x, xlab=xlab, ylab=ylab, type="n")
    points(qq, x, type="p", col=ifelse(xmiss[ix], "red", "black"))

    if(id.n > 0) {
        ind <- (n-id.n+1):n
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        text(qq[ind] + xrange/50, x[ind], ix[ind], col=ifelse(xmiss[ix], "red", "black"))
    }
    abline(0, 1, lty=2)
    title(main=eval(substitute(expression(paste(chi^2, " QQ-Plot")))))
}

##  Draw pairwise scatter plots for the data set 'x'
##  - upper triangle - scatter plot with classical and robust 0.975-ellipses
##  - histograms on the diagonal
##  - lower triangle - robust (MCD) and classical correlations
##
##  - x     - data
##  - main  - caption of the plot
##
.rrpairsNA <- function(obj, main="", sub="", xlab="", ylab="", ...){
    hcol       <- "cyan"      # colour for histogram
    dcol       <- "red"       # color of the density line
    ecol.class <- "blue"      # colour for classical ellipse
    ecol.rob   <- "red"       # colour for robust ellipse

    ## quick and dirty simulation of panel.number() which seems not to work
    ##  in the panel functions of pairs()
    ##  Retursns list of the corresponding i and j
    ##
    ##  Example call: which.ij(hbk[,1], hbk[,3], getData(CovMcd(hbk)))
    ##
    which.ij <-function(x, y, data)
    {
        i <- j <- 0
        for(k in 1:ncol(data))
        {
            ifi <- all.equal(x, data[,k], check.attributes=FALSE)
            ifj <- all.equal(y, data[,k], check.attributes=FALSE)
            if(i == 0 && !is.character(ifi) && ifi)
                i <- k
            if(j == 0 && !is.character(ifj) && ifj)
                j <- k
            if(i != 0 & j != 0)
                break
        }
        list(i=i, j=j)
    }

    panel.hist <- function(x, ...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col=hcol, ...)
    }

    panel.hist.density <- function(x,...)
    {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5) )

        h <- hist(x, plot = FALSE)
        breaks <- h$breaks; nB <- length(breaks)
        y <- h$counts; y <- y/max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col=hcol)

        tryd <- try( d <- density(x, na.rm=TRUE, bw="nrd", adjust=1.2), silent=TRUE)
        if(class(tryd) != "try-error")
        {
            d$y <- d$y/max(d$y)
            lines(d, col=dcol)
        }
    }

    panel.cor <- function(x, y, digits=2, ...)
    {
        ix <- which.ij(x, y, getData(obj))

        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))

        r <- cor(x, y, use="pairwise.complete.obs")
        rr <- getCorr(obj)[ix$i,ix$j]

        prefix <- ""
        rprefix <- ""
        ff <- 0.35

        txt  <- format(c(r, 0.123456789), digits=digits)[1]
        txt  <- paste(prefix, txt, sep="")
        rtxt <- format(c(rr, 0.123456789), digits=digits)[1]
        rtxt <- paste(rprefix, rtxt, sep="")
        txt  <- paste("(", txt, ")", sep="")

        cex  <- ff/strwidth(txt)
        text(0.5, 0.3, txt, cex = cex, col=ecol.class)
        text(0.5, 0.5, rtxt, cex = cex, col=ecol.rob)
    }
    panel.ellipse <- function(x, y, ...)
    {
        usr <- par("usr"); on.exit(par(usr))

        ix <- which.ij(x, y, getData(obj))

        cobj <- CovNAClassic(getData(obj))
        C.ls <- getCov(cobj)[c(ix$i,ix$j), c(ix$i,ix$j)]
        m.ls <- getCenter(cobj)[c(ix$i,ix$j)]
        d2.99 <- qchisq(0.975, df = 2)

        C.rr <- getCov(obj)[c(ix$i,ix$j), c(ix$i,ix$j)]
        m.rr <- getCenter(obj)[c(ix$i,ix$j)]

        e.class <- ellipsoidPoints(C.ls, d2.99, loc=m.ls)
        e.rob <- ellipsoidPoints(C.rr, d2 = d2.99, loc=m.rr)

        xmin <- min(c(min(x, na.rm=TRUE), min(e.class[,1]), min(e.rob[,1])))
        xmax <- max(c(max(x, na.rm=TRUE), max(e.class[,1]), max(e.rob[,1])))
        ymin <- min(c(min(y, na.rm=TRUE), min(e.class[,2]), min(e.rob[,2])))
        ymax <- max(c(max(y, na.rm=TRUE), max(e.class[,2]), max(e.rob[,2])))

        ff <- 0.1
        xoff <- ff*(xmax-xmin)
        yoff <- ff*(ymax-ymin)
        xmin <- xmin - xoff; #print(ff*(xmax-xmin))
        xmax <- xmax + xoff; #print(ff*(xmax-xmin))
        ymin <- ymin - yoff; #print(ff*(ymax-ymin))
        ymax <- ymax + yoff; #print(ff*(ymax-ymin))
        par(usr = c(xmin, xmax, ymin, ymax))

        acol <- ifelse(is.na(x) | is.na(y), "red", "black")
        x[is.na(x)] <- xmin
        y[is.na(y)] <- ymin

        points(x,y, col=acol)
        lines(e.class, col=ecol.class, lty="dashed")
        lines(e.rob, col=ecol.rob)
    }

    ## get the data
    x <- getData(obj)

    pairs(x, main = main, sub=sub,
        lower.panel=panel.cor,
        diag.panel=panel.hist.density,
        upper.panel=panel.ellipse,
        labels=names(getCenter(obj)),
        ...)
}

.tolellipseNA <- function(rcov, ccov, cutoff = NULL, id.n = NULL, tol = 1e-07,
    main = "Tolerance ellipse (97.5%)",
    xlab = "",
    ylab = "",
    labs,
    ...)
{
    d2.99 <- qchisq(0.975, df = 2)

    leg <- new(".Legend")
    usr <- par("usr"); on.exit(par(usr))

    if(missing(rcov) || is.null(rcov)){
        if(missing(ccov) || is.null(ccov))
            stop("Location and scatter matrix must be provided!")

        ## there is no robust location/scatter object
        rcov <- ccov
        leg@leg <- FALSE
    }
    if(is.null(data <- getData(rcov)))
        stop("No data provided!")
    n <- dim(data)[1]
    p <- dim(data)[2]
    if(p != 2)
        stop("Dimension must be 2!")

    r.cov <- getCov(rcov)
    r.loc <- getCenter(rcov)
    if(length(r.loc) == 0 ||  length(r.cov) == 0)
        stop("Invalid 'rcov' object: attributes center and cov missing!")
    z1 <- ellipsoidPoints(A=r.cov, d2=d2.99, loc=r.loc)
    rd <- sqrt(getDistance(rcov)$d)
    x1 <- c(min(data[, 1], z1[, 1], na.rm=TRUE), max(data[,1],z1[,1], na.rm=TRUE))
    y1 <- c(min(data[, 2], z1[, 2], na.rm=TRUE), max(data[,2],z1[,2], na.rm=TRUE))
    classic <- FALSE

    if(!missing(ccov) && !is.null(ccov)){
        c.cov <- getCov(ccov)
        c.loc <- getCenter(ccov)
        if(length(c.loc) == 0 ||  length(c.cov) == 0)
            stop("Invalid 'ccov' object: attributes center and cov missing!")
        classic <- TRUE
        z2 <- ellipsoidPoints(A=c.cov, d2=d2.99, loc=c.loc)

        md <- sqrt(getDistance(ccov)$d)
        x1 <- c(min(data[, 1], z1[, 1], z2[, 1], na.rm=TRUE), max(data[,1],z1[,1], z2[,1], na.rm=TRUE))
        y1 <- c(min(data[, 2], z1[, 2], z2[, 2], na.rm=TRUE), max(data[,2],z1[,2], z2[,2], na.rm=TRUE))
    }

    ## Note: the *calling* function may pass a 'missing' value
    if(missing(cutoff) || is.null(cutoff))
        cutoff <- sqrt(qchisq(0.975, df = 2))
    if(missing(id.n) || is.null(id.n))
        id.n <- sum(rd > cutoff)
    if(missing(labs) || is.null(labs))
        labs <- 1:length(rd)

    ind <- sort(rd, index.return=TRUE)$ix
    ind <- ind[(n-id.n+1):n]

    ff <- 0.1
    xoff <- ff*(x1[2]-x1[1])
    yoff <- ff*(y1[2]-y1[2])
    x1[1] <- x1[1] - xoff; #print(ff*(xmax-xmin))
    x1[2] <- x1[2] + xoff; #print(ff*(xmax-xmin))
    y1[1] <- y1[1] - yoff; #print(ff*(ymax-ymin))
    y1[2] <- y1[2] + yoff; #print(ff*(ymax-ymin))
    par(usr = c(x1[1], x1[2], y1[1], y1[2]))


##  1. Robust tolerance
##  define the plot, plot a box, plot the "good" points,
##  plot the outliers either as points or as numbers depending on outflag,
##  plot the ellipse, write a title of the plot
    plot(data, xlim = x1, ylim = y1, xlab = xlab, ylab = ylab, type="n", ...)
    box()

    ## VT:::03.08.2008
    if(id.n > 0){
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        text(data[ind, 1] + xrange/50, data[ind, 2], labs[ind])
    }

    if(leg@leg)
        points(z1, type = "l", lty=leg@lty[1], col=leg@col[1])
    title(main)

##  2. Classical tolerance ellipse and legend
    if(classic){
        points(z2, type = "l", lty=leg@lty[2], col=leg@col[2])
        if(leg@leg)
            legend("bottomright", leg@txt, pch=leg@pch, lty=leg@lty, col=leg@col)
    }

    acol <- ifelse(is.na(data[,1]) | is.na(data[,2]), "red", "black")
    data[is.na(data[,1]),1] <- x1[1]
    data[is.na(data[,2]),2] <- y1[1]

    points(data[,1], data[,2], col=acol)

    invisible()
}

## Distance plot for incomplete data:
##  Plot the robust and classical distances against the the index
##  obj - A CovRobust object,
##  getData(obj)    - data frame or matrix
##
.xydistplotNA <- function(obj, cutoff,
    main="Distance Plot",
    xlab="Index",
    ylab="Mahalanobis distance",
    col="darkred",
    colmiss="red",
    labs,
    ...)
{
    myPanel <- function(x, y, subscripts, cutoff, id.n, ...) {
        panel.xyplot(x, y, ...)
        panel.abline(h=cutoff,lty="dashed")

        n <- length(y)
        if(missing(id.n))
            id.n <- length(which(y > cutoff))
        if(id.n > 0){
            ind <- sort(y, index.return=TRUE)$ix
            ind <- ind[(n-id.n+1):n]

            xrange<-range(y)
            adj <- (xrange[2]-xrange[1])/20
            ltext(x[ind] + adj, y[ind] + adj, ind, cex=0.85)
        }
    }

    rdist <- getDistance(obj)$d
    X <- getData(obj)
    n <- nrow(X)
    p <- ncol(X)
    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))
    if(missing(labs) || is.null(labs))
        labs <- 1:length(rdist)

    acol <- ifelse(.xmiss(X), colmiss, col)

    dd1 <- sqrt(rdist)                              # robust distances
    vv<-CovNAClassic(X)
    cdist <- getDistance(vv)$d
    dd2 <- sqrt(cdist)                              # classical distances
    dd  <- c(dd1, dd2)                              # a vector with both robust and classical distances

    ind <- c(1:n, 1:n)                              # 1, 2, ..., n, 1, 2, ...n      -
    gr  <- as.factor(c(rep(1,n), rep(2,n)))         # 1, 1, ..., 1, 2, 2, ...2      -   n x 1, n x 2
    levels(gr)[1] <- "Robust"
    levels(gr)[2] <- "Classical"
    ind.col <- c(acol, acol)

    xyplot(dd~ind|gr,
                cutoff=cutoff,
                panel = myPanel,
                xlab=xlab,
                ylab=ylab,
                main=main,
                col=ind.col,
                ...)
}

## QQ-Chi plot for incomplete data:
##  Plot QQ plot of the robust and classical distances against the
##  quantiles of the Chi2 distr
##  X - data frame or matrix
##
.xyqqchi2NA <- function(obj, cutoff,
    main=eval(substitute(expression(paste(chi^2, " QQ-Plot")))),
    xlab=eval(substitute(expression(paste("Sqrt of the quantiles of the ", chi^2, " distribution")))),
    ylab="Mahalanobis distance",
    col="darkred",
    colmiss="red",
    labs,
    ...)
{
    myPanel <- function(x, y, subscripts, cutoff, id.n, ...)
    {
        y <- sort(y, index.return=TRUE)
        iy <- y$ix
        y <- y$x
        panel.xyplot(x, y, ...)
        panel.abline(0,1,lty="dashed")

        n <- length(y)
        if(missing(id.n))
            id.n <- length(which(y > cutoff))
        if(id.n > 0){
            ind <- (n-id.n+1):n

            xrange<-range(y)
            adj <- (xrange[2]-xrange[1])/50
            ltext(x[ind] + adj, y[ind] + adj, iy[ind], cex=0.85)
        }
    }

    rdist <- getDistance(obj)$d
    X <- getData(obj)
    n <- nrow(X)
    p <- ncol(X)
    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))
    if(missing(labs) || is.null(labs))
        labs <- 1:length(rdist)

    acol <- ifelse(.xmiss(X), colmiss, col)

    dd1 <- sqrt(rdist)                              # robust distances
    vv<-CovNAClassic(X)
    cdist <- getDistance(vv)$d
    dd2 <- sqrt(cdist)                              # classical distances
    dd  <- c(dd1, dd2)                              # a vector with both robust and classical distances

    qq <- sqrt(qchisq(((1:n)-1/3)/(n+1/3), p))
    ind<-c(qq, qq)
    gr<-as.factor(c(rep(1,n), rep(2,n)))            # 1, 1, ...., 1, 2, 2, ..., 2   - n x 1, n x 2
    levels(gr)[1]<-"Robust"
    levels(gr)[2]<-"Classical"
    ind.col <- c(acol, acol)

    xyplot(dd~ind|gr,
        cutoff=cutoff,
        panel = myPanel,
        xlab=xlab,
        ylab=ylab,
        main=main,
        col=ind.col,
        ...)
}
