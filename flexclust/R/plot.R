#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: plot.R 3 2013-06-12 10:06:43Z leisch $
#

setGeneric("plot")

setMethod("plot", signature(x="kcca",y="missing"),
function(x, y, which=1:2, project=NULL,
         data=NULL, points=TRUE,
         hull=TRUE, hull.args=NULL, 
         number = TRUE, simlines=TRUE,
         lwd=1, maxlwd=8*lwd, cex=1.5, numcol=FALSE, nodes=16,
         add=FALSE, xlab=NULL, ylab=NULL, xlim = NULL,
         ylim = NULL, pch=NULL, col=NULL, ...)
{
    if(length(which)!=2)
        stop(sQuote("which"), " must have length 2")

    ## try to get data from cluster object
    if(is.null(data)) data <- flexclust:::getData(x)

    ## if still NULL we cannot draw points and hulls even if user wants to
    if(is.null(data)){
        points <- hull <- FALSE
    }

    if(is.null(project))
        project <- function(x) x
    else{
        if(!is(project, "function")){
            probject <- project
            project <- function(x) Predict(probject, x)
        }
    }    
    centers <- project(x@centers)[,which]
    if(!is.null(colnames(centers))){
        if(is.null(xlab)) xlab <- colnames(centers)[which[1]]
        if(is.null(ylab)) ylab <- colnames(centers)[which[2]]
    }
    if(is.null(xlab)) xlab <- ""
    if(is.null(ylab)) ylab <- ""
    
    xrange <- apply(centers, 2, range)

    if(is(hull, "function")){
        HULLFUN <- hull
        hull <- !is.null(data)
    }
    else if(is(hull, "character")){
        hull <- match.arg(hull, c("convex", "ellipse"))
        HULLFUN <- switch(hull,
                          convex=clusterHulls,
                          ellipse=clusterEllipses)
        hull <- !is.null(data)
    }        
    else if(hull){
        HULLFUN <- clusterHulls
    }

    if(is.null(col)) col <- rep(flxColors(color="full"), length=x@k)

    if(is.logical(numcol)){
        if(numcol)
            numcol <- col
        else
            numcol <- "black"
    }
    else
        numcol <- rep(numcol, length=x@k)

    if(!is.null(data)){
        data <- x@family@preproc(data)
        distmat <- x@family@dist(data, x@centers)
        cluster <- x@family@cluster(distmat=distmat)
        cldist <- distmat[cbind(1:nrow(data), cluster)]

        data <- project(data)[,which]
        xrange <- apply(data, 2, range)
    }

    if(points){
        if(is.null(pch)){
            pch <- (cluster %% 10)
            pch[pch==0] <- 10
        }
        else if(length(pch)!=nrow(data)){
            pch <- rep(pch, length=x@k)
            pch <- pch[cluster]
        }
    }
    
    if(!add){
        if(points && !is.null(data))
            plot(data, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
                 pch=pch, col=col[cluster], ...)
        else{
            if(is.null(xlim)) xlim <- sort(xrange[,1])
            if(is.null(ylim)) ylim <- sort(xrange[,2])
            plot(0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
        }
    }
    else{
        if(points)
            points(data, pch=pch, col=col[cluster], ...)
    }

    if(hull && !is.null(data))
        do.call("HULLFUN",
                c(list(data=data,
                       cluster=cluster,
                       dist=cldist,
                       col=col),
                  hull.args))
    
    if(simlines){
        sim <- maxlwd*(x@clsim+t(x@clsim))/2

        for(k in 1:(x@k-1)){
            for(m in k:x@k)
                if(sim[k,m]>0)
                    lines(centers[c(k,m),], lwd=sim[k,m])
        }
    }

    if(number){
        rad <- ceiling(log10(x@k)) + 1.5
        points(centers, pch=21, cex=rad*cex, bg="white")
        text(centers, labels=1:x@k, cex=cex, col=numcol)
    }
    else
        points(centers, pch=nodes)

    invisible()
})

clusterEllipses <- function(data, cluster, dist, col,
                            level=c(0.5, 0.95), lwd=2*par("lwd"))
{
    level <- sort(rep(level, length=2))
    col <- rep(col, length=max(cluster))
    for(k in 1:max(cluster)){
        ok <- cluster==k
        if(sum(ok)>3){
            lines(ellipse::ellipse(cov(data[ok,]),
                                   centre=colMeans(data[ok,]),
                                   level=level[1]),
                  col=col[k], lwd=lwd)
            lines(ellipse::ellipse(cov(data[ok,]),
                                   centre=colMeans(data[ok,]),
                                   level=level[2]),
                  col=col[k], lty=2, lwd=lwd/2)
        }
    }
}

clusterHulls <- function(data, cluster, dist, col,
                         density=0, lwd=2*par("lwd"))
{
    K <- max(cluster)
    col <- rep(col, length=K)
    density <- rep(density, length=K)
    ang <- (1:K)*180/K

    cok <- complete.cases(data)
    for(k in 1:max(cluster)){
        ok <- cok & (cluster==k)
        if(length(ok)>3){
            ok1 <- ok & (dist < median(dist[ok]))
            if(sum(ok1)>3){
                hpts <- chull(data[ok1,])
                polygon(data[ok1,][c(hpts, hpts[1]),],
                        border=col[k], lwd=lwd,
                        density=density[k], col=col[k], angle=ang[k])
            }
            ok1 <- ok & (dist < 2.5*median(dist[ok]))
            if(sum(ok1)>3){
                hpts <- chull(data[ok1,])
                polygon(data[ok1,][c(hpts, hpts[1]),], border=col[k],
                        lty=2, lwd=lwd/2)
            }
        }
    }
}



###**********************************************************

setGeneric("image")

setMethod("image", signature(x="kcca"),
function(x, which = 1:2, npoints = 100,
         xlab = "", ylab = "", fastcol = TRUE, col=NULL,
         clwd=0, graph=TRUE, ...) 
{
    if(length(which)!=2)
        stop(sQuote(which), "must be a vector of length 2")
    
    X <- seq(x@xrange[1,which[1]],
             x@xrange[2,which[1]], length=npoints)
    Y <- seq(x@xrange[1,which[2]],
             x@xrange[2,which[2]], length=npoints)
    Y1 <- rep(Y, rep.int(length(X), length(Y)))
    X1 <- rep(X, times = ceiling(length(Y1)/length(X)))
    
    neib1 <- neighbours(x)

    if(fastcol){
        Z <- x@family@cluster(cbind(X1,Y1), x@centers[,which], n=1)
        neib2 <- neib1
    }
    else{
        Z <- x@family@cluster(cbind(X1,Y1), x@centers[,which], n=2)
        neib2 <- neighbours(Z)
        Z <- Z[[1]]
    }

    Z <- matrix(Z, nrow=length(X))
    image(X, Y, Z,
          col=neighbourColors(neib2, col=col), xlab=xlab, ylab=ylab)

    if(clwd>0)
        contour(X, Y, Z, levels=seq(1.5, max(Z)),
                add=TRUE, drawlabels=FALSE, lwd=clwd)

    if(graph)
        plot(x, which=which, add=TRUE, ...)
    
    invisible(list(x=X,y=Y,z=Z))
})



neighbourColors <- function(object, col=NULL)
{
    if(is.null(col))
        col <- LightColors
    
    icol <- seq(along=col)
    z <- rep(NA, length(object))
    names(z) <- 1:length(z)

    nmis <- 0
    while(any(is.na(z))){
        
        node <- min(which(is.na(z)))
        nodes <- c(node,object[[node]])
        
        for(n in nodes){
            if(is.na(z[n])){
                used <- na.omit(z[object[[n]]])
                if(length(used)>=length(col)){
                    nmis <- nmis+1
                    z[n] <- sample(icol, 1)
                }
                else if(length(used>0))
                    z[n] <- icol[-used][1]
                else
                    z[n] <- icol[1]
            }
        }
    }
    if(nmis>0)
        warning("colors were missing ", nmis, " times\n")
    
    col[z]
}

###**********************************************************

setMethod("plot", signature(x="stepFlexclust", y="missing"),
function(x, y, type=c("barplot", "lines"), totaldist=NULL,
          xlab=NULL, ylab=NULL, ...)
{
    type <- match.arg(type)
    X <- x@k
    Y <- sapply(x@models, function(z) info(z, "distsum"))

    if(is.null(totaldist))
        totaldist <- 2 %in% X

    if(totaldist){
        X <- c(1, X)
        Y <- c(x@totaldist, Y)
    }
        
    if(is.null(xlab))
        xlab <- "number of clusters"
    if(is.null(ylab))
         ylab <- "sum of within cluster distances"

    if (type == "barplot") 
        barplot(Y, names = X,
                xlab=xlab, ylab = ylab, ...)
    else {
        plot(X, Y, type = "b", axes = FALSE,
             xlab = xlab, ylab = ylab, ...)
        axis(2)
        axis(1, at = X, labels = X)
    }

    invisible()
})

###**********************************************************

setGeneric("pairs")

setMethod("pairs", signature(x="kcca"),
function(x, which=NULL, project=NULL, oma=NULL, ...)
{
    if(is.null(which))
        which <- 1:ncol(x@centers)

    if(is.null(oma))
        oma <- rep(3,4)

    if(is.null(project)){
        NAMES <- colnames(x@centers)[which]
    }
    else{
        pc <- Predict(project, x@centers)
        NAMES <- colnames(pc)[which]
    }

    n <- length(which)

    opar <- par("mar", "oma", "mfrow")
    on.exit(par(opar))

    par(mfcol=c(n, n), oma=oma, mar=rep(0.2,4))

    for(k in 1:n){
        for(l in 1:n){

            if(k!=l){
                plot(x, which=which[c(k,l)], project=project, axes=FALSE, ...)
                if(l==1) axis(side=3)
                if(l==n) axis(side=1)
                if(k==1) axis(side=2)
                if(k==n) axis(side=4)
            }
            else{
                plot.new()
                plot.window(0:1, 0:1)
                text(0.5, 0.5, NAMES[k], cex=2)
            }
            box()
        }
    }
    
})

###**********************************************************

setMethod("plot", signature(x="kccasimple",y="missing"),
function(x, y)
{
    stop("No plot method for kccasimple objects available!\n")
})

setMethod("image", signature(x="kccasimple"),
function(x)
{
    stop("No image method for kccasimple objects available!\n")
})

## setMethod("barplot", signature(height="kccasimple"),
## function(height)
## {
##     warning("No barplot method for kccasimple objects available!\n")
## })

setMethod("pairs", signature(x="kccasimple"),
function(x)
{
    stop("No pairs method for kccasimple objects available!\n")
})

###**********************************************************

### A wrapper for predict to allow special cases like lda
Predict <- function(object, newdata)
{
    if(inherits(object, "lda"))
        return(predict(object, newdata)$x)
    else
        return(predict(object, newdata))
}
