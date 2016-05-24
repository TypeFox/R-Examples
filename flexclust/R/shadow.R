#
#  Copyright (C) 2005 Friedrich Leisch
#  $Id: shadow.R 3 2013-06-12 10:06:43Z leisch $
#

setClass("shadow",
         representation(values="numeric",
                        cluster="integer",
                        size="integer",
                        k="integer",
                        similarity="logical"))

setGeneric("shadow", function(object, ...) standardGeneric("shadow"))

setMethod("shadow", signature(object="kccasimple"),
function(object, ...)
{
    z <- (2*object@cldist[,1])/rowSums(object@cldist[,1:2])

    new("shadow", values=z, cluster=object@cluster, k=object@k,
        size=info(object, "size"), similarity=TRUE)
})

setMethod("show", "shadow",
function(object)
{
    z <- tapply(object@values, object@cluster, mean)
    print(z)
    invisible(z)
})

setMethod("plot", "shadow", function(x, y, ...) shadowPlot(x, ...))

###**********************************************************

setGeneric("shadowPlot", function(object, ...)
    standardGeneric("shadowPlot"))

setMethod("shadowPlot", signature(object="shadow"),
function(object, nrow=NULL, varwidth=is.null(nrow),
         panel=panelShadow, ylim=NULL, ...)
{
    y <- object@values

    if(is.null(ylim))
        ylim <- c(min(0, min(y)), 1)

    grid.newpage()
    pushViewport(viewport(width=0.95, height=0.95))

    if(varwidth){
        nrow <- 1
        ncol <- object@k
        pushViewport(viewport(layout=grid.layout(nrow=nrow, ncol=ncol,
                              widths = object@size/sum(object@size))))
    }
    else{
        if(is.null(nrow)){
            ncol <- ceiling(sqrt(object@k))
            nrow <- ceiling(object@k / ncol)
        }
        else{
            ncol <- ceiling(object@k / nrow)
        }
        pushViewport(viewport(layout=grid.layout(nrow, ncol)))
    }

    k <- 1
    for(i in 1:nrow) for(j in 1:ncol) {
        if(k > object@k) invisible(y)
        
        yk <- sort(y[object@cluster==k], decreasing=TRUE)
        
        pushViewport(viewport(layout.pos.row = i,
                              layout.pos.col = j))
        grid.rect()

        pushViewport(viewport(x=0, y=1, height=unit(1, "lines"),
                              just=c("left", "top")))
        grid.rect()
        grid.text(k)
        popViewport()
        
        pushViewport(viewport(x=0, y=0,
                              height=unit(1, "npc")-unit(1, "lines"),
                              just=c("left", "bottom")))
        panel(yk, ylim=ylim, breaks=min(length(yk), 100),
              similarity=object@similarity, ...)
        popViewport(2)
        k <- k+1
    }
    popViewport()
})

setMethod("shadowPlot", signature(object="kccasimple"),
function(object, ...) shadowPlot(shadow(object), ...))

###**********************************************************

setGeneric("Silhouette", function(object, ...)
    standardGeneric("Silhouette"))

setMethod("Silhouette", signature(object="kcca"),
function(object, data=NULL, ...)
{
    if(is.null(data)){
        data <- getData(object)
        cluster <- object@cluster
    }
    else {
        data <- as.matrix(data)
        cluster <- object@family@cluster(data, object@centers)
    }
    
    z <- rep(1, length(cluster))

    for(k in 1:object@k)
    {
        ok1 <- cluster==k
        if(!any(ok1)) next

        a <- rowMeans(object@family@dist(data[ok1,,drop=FALSE],
                                         data[ok1,,drop=FALSE]))
        b <- rep(Inf, length(a))
        
        for(m in (1:object@k)[-k]){
            ok2 <- cluster==m
            if(!any(ok2)) next
            
            b <- pmin(b, rowMeans(object@family@dist(data[ok1,,drop=FALSE],
                                                     data[ok2,,drop=FALSE])))

            z[ok1] <- (b-a) / pmax(a,b)
        }
    }
                                    
    new("shadow", values=z, cluster=cluster, k=object@k,
        size=as.vector(table(cluster)), similarity=FALSE)
})



###**********************************************************

panelShadow <- function(x, y=NULL, breaks=length(x),
                        xlim=NULL, ylim=0:1,
                        col1=NULL, col2=NULL,
                        similarity=TRUE)
{
    if(is.null(y)){
        y <- x
        x <- seq(0,1,along=y)
    }

    if(is.null(col1))
        col1 <- hcl(10, c=c(50, 30), l=c(55, 85))
    
    if(is.null(col2))
        col2 <- hcl(130, c=c(50, 30), l=c(55, 85))

    if(!similarity){
        col3 <- col1
        col1 <- col2
        col2 <- col3
    }
    
    ym <- mean(y)

    ## approx garantiert sortierte x-werte -> auch benutzen wenn nicht
    ## notwendig, d.h. breaks=length(x)
    y <- approx(x, y, n=breaks)
    x <- c(min(x), y$x, max(x))
    y <- c(0, y$y, 0)
    
    if(is.null(xlim)) xlim <- range(x)
    
    pushViewport(viewport(xscale=xlim, yscale=ylim))
    
    grid.segments(x0=unit(0, "npc"),
                  y0=unit(c(-.5, .5), "native"),
                  x1=unit(1, "npc"),
                  y1=unit(c(-.5, .5), "native"),
                  gp=gpar(col="grey", lty=2))
                  
    grid.segments(x0=unit(0, "npc"),
                  y0=unit(seq(-.75,.75,by=0.5), "native"),
                  x1=unit(1, "npc"),
                  y1=unit(seq(-.75,.75,by=0.5), "native"),
                  gp=gpar(col="grey", lty=3))
                  
    if(ym>0)
        grid.rect(x=0, y=0, width=1, height=ym, default.units="native",
                  just=c("left","bottom"), gp=gpar(fill=col1[2]))
    else
        grid.rect(x=0, y=0, width=1, height=ym, default.units="native",
                  just=c("left","bottom"), gp=gpar(fill=col2[2]))

    if(ylim[1]>=0){
        grid.polygon(unit(x, "native"),
                     unit(y, "native"),
                     gp=gpar(fill=col1[1]))
    }
    else if(ylim[2]<=0){
        grid.polygon(unit(x, "native"),
                     unit(y, "native"),
                     gp=gpar(fill=col2[1]))
    }
    else{
        pushViewport(viewport(x=0.5, y=1,
                              xscale=xlim,
                              yscale=c(0, ylim[2]),
                              height=unit(ylim[2], "native"),
                              just=c("centre", "top"),
                              clip="on"))
        grid.polygon(unit(x, "native"),
                     unit(y, "native"),
                     gp=gpar(fill=col1))
        popViewport()
        
        pushViewport(viewport(x=0.5, y=0,
                              xscale=xlim,
                              yscale=c(ylim[1], 0),
                              height=unit(-ylim[1], "native"),
                              just=c("centre", "bottom"),
                              clip="on"))
        grid.polygon(unit(x, "native"),
                     unit(y, "native"),
                     gp=gpar(fill=col2))
        popViewport()
    }

##     grid.segments(x1=unit(1,"npc")-unit(rep(c(1, 0.5), length=9), "char"),
##                   y0=unit(seq(-1,1,by=0.25), "native"),
##                   x0=unit(1, "npc"),
##                   y1=unit(seq(-1,1,by=0.25), "native"))
                  

    popViewport()
}

    
###**********************************************************

shadowStars <- function(object, which=1:2, project=NULL,
                        width=1, varwidth=FALSE,
                        panel=panelShadowStripes, box=NULL,
                        col=NULL, add=FALSE, ...)
{
    if(is.null(project))
        centers <- object@centers
    else{
        if(is(project, "function"))
            centers <- project(object@centers)
        else
            centers <- predict(project, object@centers)
    }

    col <- expandColors(col, object)
    
    if(!add) grid.newpage()

    pushViewport(viewport(width=unit(0.9, "snpc"),
                          height=unit(0.9, "snpc"),
                          xscale=range(centers[,which]),
                          yscale=range(centers[,which]),
                          name="FL"))

    N <- length(object@cluster)
    d <- 2*object@cldist[,1]/rowSums(object@cldist)

    width <- width*max(dist(centers))/25
    maxsize <- max(info(object, "size"))

    for(k in 1:object@k){
        for(m in (1:object@k)[-k]){
            if(any(ok <- (object@cluster==k & object@second==m))){

                x <- centers[c(k,m),which]
                breite <- dist(x)[1]/2
                hoehe <- width
                if(varwidth)
                    hoehe <-  2 * hoehe * sum(ok)/maxsize
                
                winkel <- 180*atan2(diff(x[,2]),diff(x[,1]))/pi
                
                pushViewport(viewport(x=x[1,1], y=x[1,2],
                                      height=hoehe,
                                      width=breite, angle=winkel,
                                      default.units="native",
                                      just=c("left", "center"),
                                      yscale=c(-1,1)))

                if(!is.null(box))
                    grid.rect(gp=gpar(col=box))
                panel(d[ok], col=col[ok], ...)
                popViewport()
            }
        }
    }

    rad <-
        if(object@k<10) unit(0.8, "char")
        else unit(1, "char")
    
    grid.circle(x=centers[,which[1]],
                y=centers[,which[2]],
                r=rad,
                default.units="native", gp=gpar(fill="white"))
        
    grid.text(x=centers[,which[1]],
              y=centers[,which[2]],
              label=1:object@k,
              default.units="native", gp=gpar(col="red"))
    
    popViewport()
}

panelShadowViolin <- function(x, ...)
{
    if(length(x)<=4) return()
    
    d <- density(x, n=100, from=0, to=1)
    x <- d$x
    ## divide by 2 seems to give nice default height
    y <- d$y/2
    
    grid.polygon(c(x, rev(x)), unit(c(y, -rev(y)), "native"),
                 gp=gpar(fill="grey70"), name="FL")
}

panelShadowBP <- function(x, ...)
{
    if(length(x)<=4) return()
    
    y <- rev(seq(0,1, 0.025))
    x <- c(quantile(x, y), 0)
    y[y>0.5] <- 1 - y[y>0.5]
    ## multiply by 2 gives similar width to density in violin
    y <- 2*c(y, 0)

    grid.polygon(c(x, rev(x)), unit(c(y, -rev(y)), "native"),
                 gp=gpar(fill="grey90"), name="FL")

    xp <- x[c("25%", "50%", "75%")]
    grid.segments(xp, unit(c(.5, 1, .5), "native"),
                  xp, unit(-c(.5, 1, .5), "native"))

}

panelShadowSkeleton <- function(x, ...){

    m <- mean(x)
    grid.rect(height=m, gp=gpar(fill="grey70"))

}

panelShadowStripes <- function(x, col, ...)
{
    grid.segments(x0=x, x1=x, y0=0, y1=1, gp=gpar(col=col))
}

###**********************************************************

