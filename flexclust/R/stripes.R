stripes <- function(object, groups=NULL, type=c("first", "second", "all"),
                    beside=(type!="first"), col=NULL, gp.line=NULL, gp.bar=NULL,
                    gp.bar2=NULL, number=TRUE, legend=!is.null(groups),
                    ylim=NULL, ylab="distance from centroid",
                    margins=c(2,5,3,2), ...)
{
    type <- match.arg(type)
    k <- object@k

    xd <- NULL
    if(type=="all"){
        xd <- object@family@dist(flexclust:::getData(object),
                                 object@centers)
    }
    
    if(is.null(ylim)){
        ylim <- c(0, maxDist(object, type, xd))
    }
    yticks <- pretty(ylim)
    ym <- max(yticks)
    yr <- range(yticks)

    if(!is.null(groups)){
        groups <- as.factor(groups)        
        ng <- length(levels(groups))
        if(is.null(col)) col <- 1:ng
        col <- rep(col, length=ng)

        if(legend){
            legend <- list(labels=levels(groups), col=col, pch=15)
        }
        else{
            legend <- NULL
        }
        
        col <- col[as.integer(groups)]            
    }
    else{
        col <- expandColors(col, object)
        legend <- NULL
    }

    if(is.null(gp.line)) gp.line <- gpar()

    if(is.null(gp.bar)) gp.bar=gpar(col=grey(.6), fill=grey(.75))
    if(is.null(gp.bar2)) gp.bar2=gpar(col=grey(.8), fill=grey(.95))

    pf <- plotFramework(ylim=yr, margins=margins, xaxis=FALSE, yaxis=TRUE,
                        ylab=ylab, legend=legend, ...)
    
    pushViewport(viewport(layout=grid.layout(nrow=1, ncol=k)))

    m <- info(object, "max_dist")
    
    for(n in 1:k){
        pushViewport(viewport(layout.pos.row = 1,
                              layout.pos.col = n))
        
        pushViewport(viewport(1, 0, width=0.90,
                              just=c(1,0), yscale=yr))

        if(number) grid.text(n, y=unit(-0.75, "lines"))

        x0 <- 0; w <- 1
        if(beside){
            x0 <- (n-1)/k
            w <- 1/k
            grid.rect(0, 0, height=unit(m[n], "native"),
                      just=c(0,0), gp=gp.bar2)
        }
            
        grid.rect(x0, 0, height=unit(m[n], "native"), width=w,
                  just=c(0,0), gp=gp.bar)
        if(type != "first"){
            grid.rect(0, unit(m[n], "native"),
                      height=unit(1, "npc")-unit(m[n], "native"),
                      just=c(0,0), gp=gp.bar2)
        }

        z <- getDistColClust(object, n, type, col, xd)
        
        x1 <- 1
        if(beside){
            x0 <- (z$cluster-1)/k
            x1 <- z$cluster/k
        }

        gp.line$col <- z$col

        grid.segments(x0=x0, x1=x1, y0=z$dist, y1=z$dist, gp=gp.line,
                      default.units="native")


        popViewport(2)
    }
    popViewport(2)
}


getDistColClust <- function(object, k, type, col, xd=NULL)
{
    if(type=="all"){
        z <- list(dist=xd[,k], col=col, cluster=object@cluster)
    }
    else{
        ok <- object@cluster==k
        x <- object@cldist[ok,1]
        xcol <- col[ok]
        clus <- object@cluster[ok]

        if(type=="second"){
            ok <- object@second==k
            if(any(ok)){
                x <- c(x, object@cldist[ok,2])
                xcol <- c(xcol, col[ok])
                clus <- c(clus, object@cluster[ok])
            }
        }
        z <- list(dist=x, col=xcol, cluster=clus)
    }
    z
}



maxDist <- function(object, type, xd=NULL)
{
    if(type=="first"){
        z <- max(object@cldist[,1])
    }
    else if(type=="second"){
        z <- max(object@cldist)
    }
    else{
        if(is.null(xd))
            stop("Need xd: matrix of cluster distances!\n")
        z <- max(xd)
    }
    z
}

plotFramework <- function(xlim=0:1, ylim=0:1,
                          margins=c(4.1, 5.1, 2.1, 2.1),
                          xaxis=TRUE, yaxis=TRUE,
                          xticks=NULL, yticks=NULL,
                          xlab="", ylab="", main="",
                          add=FALSE, legend=NULL)
{
    if(!add) grid.newpage()

    grid.text(main, y=unit(1, "npc") - unit(0.75, "lines"),
              gp=gpar(cex=1.5, fontface="bold"))

    if(!is.null(legend)){

        legend$col <- rep(legend$col, length=length(legend$labels))
        legend$pch <- rep(legend$pch, length=length(legend$labels))

        w <- unit( max(6, max(nchar(legend$labels)))+2, "char")
        margins[4] <- convertUnit(w, "lines")

        lx <- unit(1, "npc") - unit(margins[4]-2, "lines")
        for (n in 1:length(legend$labels)) {
            ly <- unit(1, "npc") - unit(margins[3]+n, "lines")
            grid.points(lx, ly,
                        pch=legend$pch[n],
                        gp=gpar(col=legend$col[n]))
            
            grid.text(legend$labels[n], lx+unit(1, "char"), ly, just=0)
        }
    }
    
    vp <- pushViewport(plotViewport(margins=margins, yscale=ylim))

    if(xaxis) grid.xaxis(xticks)
    if(yaxis) grid.yaxis(yticks)

    grid.text(xlab, y=unit(-3, "lines"))
    grid.text(ylab, x=unit(-4, "lines"), rot=90)

    invisible(vp)
}
