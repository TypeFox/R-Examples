setClass("projAxes", representation(arrows="list", text="list"))

projAxes <- function(object, which=1:2, center=NULL,
                     col="red", radius=NULL,
                     minradius=0.1, textargs=list(col=col),
                     col.names=getColnames(object),
                     which.names="", group=NULL, groupFun=colMeans,
                     plot=TRUE, ...)
{
    pu = par("usr")
    if(is.null(center)){
        center=c((pu[1]+pu[2])/2, (pu[3]+pu[4])/2)
    }
    else
        center = rep(center, length=2)

    D <- rbind(0, diag(length(col.names)))
    colnames(D) <- col.names
    z <- Predict(object, D)[,which]
    z0 <- z[1,]
    z <- z[-1,]
    z[,1] <- z[,1]-z0[1]
    z[,2] <- z[,2]-z0[2]

    if(!is.null(group)){
        group <- rep(group, length=nrow(z))
        gu <- unique(group)
        gu <- gu[!is.na(gu)]
        z1 <- matrix(0, nrow=length(gu), ncol=2)
        for(g in 1:length(gu)){
            ok <- group==gu[g]
            ok[is.na(ok)] <- FALSE
            z1[g,] <- groupFun(z[ok,,drop=FALSE])
        }
        z <- z1
        col.names <- as.character(gu)
    }
    
    if(is.null(radius)){
        radius <- c(abs(center[1]-pu[1:2]), abs(center[2]-pu[3:4]))
        radius <- 0.75*min(radius/apply(z,2,function(x) max(abs(x))))
    }
    
    x0 = center[1]
    y0 = center[2]
    x1 = z[,1]*radius
    y1 = z[,2]*radius

    r <-  sqrt(x1^2+y1^2)
    ok <- ( r >= (minradius * max(r)) ) &
          ( regexpr(which.names, col.names) > 0)

    x1=x1[ok]
    y1=y1[ok]

    alpha <- atan2(y1, x1)
    pos <- rep(4, length(alpha))
    pos[(alpha >= pi/4) & (alpha <= 3*pi/4)] <- 3
    pos[abs(alpha)>3*pi/4] <- 2
    pos[(alpha <= -pi/4) & (alpha>= -3*pi/4)] <- 1
    
    z <- new("projAxes",
             arrows=c(list(x0=x0, y0=y0, x1=x1+x0, y1=y1+y0, col=col), ...),
             text=c(list(x=x1+x0, y=y1+y0, labels=col.names[ok],
             pos=pos, offset=rep(0.5, length(pos))),
             textargs))
    if(plot) plot(z)
    invisible(z)
}

setMethod("plot", signature(x="projAxes", y="missing"),
function(x, y, ...)
{          
    do.call("arrows", x@arrows)
    do.call("text", x@text)
})


setGeneric("placeLabels",
function(object) standardGeneric("placeLabels")) 

setMethod("placeLabels", signature(object="projAxes"),
function(object)
{
    plot(object)
    for(n in 1:length(object@text$labels)){
        points(object@arrows$x1[n], object@arrows$y1[n], col="blue",
               cex=3, lwd=3)
        x <- locator(type="p", col="blue", pch=3, cex=2, lwd=2)
        points(object@arrows$x1[n], object@arrows$y1[n], col="white",
               cex=3, lwd=3)
        if(length(x$x)){
            points(x$x, x$y, type="p", col="white", pch=3, cex=2,
                   lwd=2)
            object@text$x[n] <- tail(x$x, n=1)
            object@text$y[n] <- tail(x$y, n=1)
            object@text$offset[n] <- 0
            object@text$pos[n] <-
                ifelse(object@text$x[n] > object@arrows$x1[n], 4, 2)
        }
        text(x=object@text$x[n], y=object@text$y[n],
             labels=object@text$labels[n], pos=object@text$pos[n],
             col="blue", offset=0)
2    }
    object
})
        

getColnames <- function(object) UseMethod("getColnames")

getColnames.prcomp <- function(object) rownames(object$rotation)

getColnames.lda <- function(object) colnames(object$means)


