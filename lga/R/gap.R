gap <- function(x, K, B, ...) UseMethod("gap")

"gap.default" <- function(x, K, B, criteria=c("tibshirani", "DandF","none"), nnode=NULL, scale=TRUE, ...){
    ## Setup the criteria
    criteria <- match.arg(criteria)
    doall <- ifelse(criteria=="tibshirani", FALSE, TRUE)

    ## Scale the dataset (if required), otherwise parse to a matrix
    x <- as.matrix(x)
    if (any(is.na(x))) stop("Missing data in x")
    if(scale)
        x <- scale(x, center=FALSE, scale=sqrt(apply(x,2,var)))
    if (!is.numeric(x)) stop("Data does not appear to be numeric.\n")
    n <- nrow(x); d <- ncol(x)

    GAP <- logWks <- rep(NA, K)
    ElogWks <- matrix(NA, nrow=K, ncol=2,
                      dimnames=list(NULL, c("ElogWks", "Sks")))
    finished <- FALSE

    ## Check that K is appropriate
    if (K > floor(n/d)) {
        cat("Not enough observations to consider",K, "clusters.Using K =",
            floor(n/d),"instead.\n" )
        K <- floor(n/d)
    }
    k <- 0

    if(!is.null(nnode)) {
        ## set up the nodes
        cat("Setting up nodes \n")
        if (!(require(snow))) stop("Can't find required packages: snow")
        cl <- makeCluster(nnode)
        clusterEvalQ(cl, library(lga))
        ## Randomize the seed (dirty, yes, but doesn't require SPRNG)
        clusterApply(cl, runif(length(cl), max=10e6),set.seed)
    }

    ## start looping through k
    while (!finished & k < K){
        k <- k+1

        ## Calculate the GAP for actual data
        cat("\nCalculating GAP at k =", k,"\nCalculating log(W_k)")
        logWks[k] <- gap.logW(x, k)

        ## Now bootstrap on reference distribution
        cat("\nCalculating E log(W_k) (Bootstrap)")

        if (is.null(nnode))
            BootOutput <- gap.boot(x, B, n, k)
        else {
            B <- ceiling(B/nnode)*nnode
            BootOutput <- unlist(clusterCall(cl, gap.boot, x, B/nnode, n, k))
        }

        ## Calculate GAP statistic, and return results
        ElogWks[k,] <- c(mean(BootOutput), sqrt(var(BootOutput)*(1+1/B)))
        GAP[k] <- ElogWks[k,1] - logWks[k]
        if (k > 1)
            if(GAP[k-1] >= GAP[k]-ElogWks[k,2] & !doall)
                finished <- TRUE
    }

    if (!is.null(nnode)){
        cat("Closing nodes\n")
        stopCluster(cl)
    }

    outdata <- cbind(Gap=GAP[1:k], logWks=logWks[1:k], ElogWks[1:k,])
    rownames(outdata) <- paste("k=", 1:dim(outdata)[1], sep="")
    output <- list(finished=finished, nclust=k-1, data=outdata, criteria=criteria)
    class(output) <- "gap"
    output$nclust <- criteria(output)


    plot(output)
    print(output)
    invisible(output)
}


plot.gap <- function(x, ...){
    GapData <- x$data
    logWk <- GapData[,2]
    ElogWk <- GapData[,3:4]
    op <- par(mfrow=c(1,2)); on.exit(par(op))

    ## First plot  - number of clusters vs Obs and Exp log(Wk)
    ysize <- range(logWk, ElogWk[,1])
    ysize <- c(floor(ysize[1]), ceiling(ysize[2]))
    yseq <- seq(from=ysize[1], to=ysize[2], length=5)
    plot(logWk, type="b", pch="O", ylim=ysize, axes=F,
         xlab="Number of clusters k", ylab="Obs and Exp log(Wk)")
    points(ElogWk[,1], type="b", pch="E")
    box()
    axis(1, 1:length(logWk), at=1:length(logWk))
    axis(2, yseq, at=yseq)

    ## Second plot - number of clusters vs Gap statistic
    plot(GapData[,1], type="b", axes=F, xlab="Number of clusters k", ylab="Gap")
    axis(1, 1:length(logWk), at=1:length(logWk))
    ysize <- range(GapData[,1], GapData[,1]-ElogWk[,2])
    ysize <- c(floor(ysize[1]), ceiling(ysize[2]))
    yseq <- seq(from=ysize[1], to=ysize[2], length=5)
    axis(2, yseq, at=yseq)
    arrows(1:dim(GapData)[1], GapData[,1], 1:dim(GapData)[1], GapData[,1] - ElogWk[,2],
           angle=90, code=2 , length=0.1)
    if(!is.na(x$nclust))
        points(x$nclust, GapData[x$nclust,1], pch=19, col='red')
}

print.gap <- function(x, ...){
    cat("\n\nCriteria used:", x$criteria,"\n")
    if (!x$criteria == "none"){
        if (is.na(x$nclust))
            cat("Nothing conclusive found for K =", nrow(x$data), "- consider raising K (if possible).\n")
        else
            cat("Gap suggests there are", x$nclust,"clusters\n\n")
    }
    print(x$data)
}

"gap.boot" <- function(xsc, B, n, k){
    ## this function generates B samples of size n using the function gap.genBox
    ## these are then given to gap.logW.
    ## Uses the boot function from the library of the same name.
    ## This returns a matrix (t) with the statistic for each replicate
    return(boot(data=xsc, statistic= gap.logW, R=B, sim="parametric",
                ran.gen=gap.genBox, mle=n, k=k)$t)
}

"gap.genBox" <- function(x, n){
    ## generates a uniform box over the range of x, transformed to the
    ## eigenvectors of x.
    y <- scale(x, scale=FALSE)
    svdOut <- svd(y)
    Ranges <- apply(y%*%svdOut$v,2,range)
    z <- apply(Ranges, 2, function(x, nn)
               runif(nn, min=x[1], max=x[2]), nn=n)
    zPrime <- z %*% t(svdOut$v)
    return(sweep(zPrime, 2, attr(y, 'scaled:center'), FUN="+"))
}


"gap.logW" <- function(x, k){
    cat(".")
    n <- dim(x)[1]
    d <- dim(x)[2]

    if (k==1)
        groups <- rep(1, n)
    else
        groups <- lga(x, k, niter=20, scale=FALSE, silent=TRUE)$cluster

    hpcoef <- matrix(NA, nrow=k, ncol=d+1)
    for (i in 1:k)
        hpcoef[i,] <- lga.orthreg(x[groups==i,])
    return(log(lga.calculateROSS(hpcoef, x, n, d, groups)))
}

"criteria" <- function(x){
    switch(x$criteria,
           none = criteria.none(x),
           DandF = criteria.DandF(x),
           tibshirani = criteria.tibshirani(x))
}

"criteria.none" <- function(x) {
    ## just return the GAP values
    return(NA)
}

"criteria.DandF" <- function(x){
    ## the method given in Dudoit and Fridlyand (2002)
    y <- x$data
    crit <- diff(y[which.max(y[,"Gap"]), c("Sks", "Gap")])
    nclust <- min(which(y[,"Gap"] > crit))
    return(ifelse(nclust == nrow(y), NA, nclust))
}

"criteria.tibshirani" <- function(x) {
    ## the method given in Tibshirani, Walter & Hastie (2001)
    return(ifelse(!x$finished, NA, x$nclust))
}
