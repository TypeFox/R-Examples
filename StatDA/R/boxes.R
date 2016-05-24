boxes <-
function(x, xA=1,yA=2,zA=3,
    labels = dimnames(x)[[1]], locations = NULL, nrow = NULL, ncol = NULL, 
    key.loc = NULL, key.labels = dimnames(x)[[2]], key.xpd = TRUE, xlim = NULL,
    ylim = NULL, flip.labels = NULL, len=1, leglen=1,
    axes = FALSE, frame.plot = axes, main = NULL,
    sub = NULL, xlab = "", ylab = "", cex = 0.8, lwd = 0.25,
    lty = par("lty"), xpd = FALSE, mar = pmin(par("mar"), 1.1 +
        c(2 * axes + (xlab != ""), 2 * axes + (ylab != ""), 1,
            0)), add = FALSE, plot = TRUE, ...)
{
# draw boxes as multivariate graphics
#
# x ... multivariate data in form of matrix or data frame
# xA, yA, zA ... assignment of clusters to the coordinates of the boxes
# labels ... vector of character strings for labeling the plots
# locations ... locations for the boxes on the plot (e.g. X/Y coordinates)
# nrow, ncol ... integers giving the number of rows and columns to use when
#          'locations' is 'NULL'.  By default, 'nrow == ncol', a square
#          layout will be used.
# key.loc ... vector with x and y coordinates of the unit key.
# key.labels: vector of character strings for labeling the segments of
#         the unit key.  If omitted, the second component of
#         'dimnames(x)' is used, if available.
# key.xpd: clipping switch for the unit key (drawing and labeling), see
#          'par("xpd")'.
# xlim: vector with the range of x coordinates to plot.
# ylim: vector with the range of y coordinates to plot.
# flip.labels: logical indicating if the label locations should flip up
#          and down from diagram to diagram. Defaults to a somewhat
#          smart heuristic.
# len, leglen: multiplicative values for the space of the labels on the legend
# axes: logical flag: if 'TRUE' axes are added to the plot.
# frame.plot: logical flag: if 'TRUE', the plot region is framed.
#   main: a main title for the plot.
#    sub: a sub title for the plot.
#   xlab: a label for the x axis.
#   ylab: a label for the y axis.
#    cex: character expansion factor for the labels.
#    lwd: line width used for drawing.
#    lty: line type used for drawing.
#    xpd: logical or NA indicating if clipping should be done, see
#         'par(xpd = .)'.
#    mar: argument to 'par(mar = *)', typically choosing smaller
#         margings than by default.
#    add: logical, if 'TRUE' _add_ boxes to current plot.
#   plot: logical, if 'FALSE', nothing is plotted.
#    ...: further arguments, passed to the first call of 'plot()', see
#         'plot.default' and to 'box()' if 'frame.plot' is true.



    if (is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        stop("'x' must be a matrix or a data frame")
    if (!is.numeric(x))
        stop("data in 'x' must be numeric")
    n.loc <- nrow(x)
    n.seg <- ncol(x)
    if (is.null(locations)) {
        if (is.null(nrow))
            nrow <- ceiling(if (!is.numeric(ncol)) sqrt(n.loc) else n.loc/ncol)
        if (is.null(ncol))
            ncol <- ceiling(n.loc/nrow)
        if (nrow * ncol < n.loc)
            stop("nrow * ncol <  number of observations")
        ff <- if (!is.null(labels)) 2.3
        else 2.1
        locations <- expand.grid(ff * 1:ncol, ff * nrow:1)[1:n.loc, ]
        if (!is.null(labels) && (missing(flip.labels) || !is.logical(flip.labels)))
            flip.labels <- ncol * mean(nchar(labels, type = "c")) > 30
    }
    else {
        if (is.numeric(locations) && length(locations) == 2) {
            locations <- cbind(rep.int(locations[1], n.loc),
                rep.int(locations[2], n.loc))
            if (!missing(labels) && n.loc > 1)
                warning("labels do not make sense for a single location")
            else labels <- NULL
        }
        else {
            if (is.data.frame(locations))
                locations <- data.matrix(locations)
            if (!is.matrix(locations) || ncol(locations) != 2)
                stop("'locations' must be a 2-column matrix.")
            if (n.loc != nrow(locations))
                stop("number of rows of 'locations' and 'x' must be equal.")
        }
        if (missing(flip.labels) || !is.logical(flip.labels))
            flip.labels <- FALSE
    }
    xloc <- locations[, 1]
    yloc <- locations[, 2]

    x[is.na(x)] <- 0
    mx <- max(x <- x * len)
    if (is.null(xlim))
        xlim <- range(xloc) + c(-mx, mx)
    if (is.null(ylim))
        ylim <- range(yloc) + c(-mx, mx)
    op <- par(mar = mar, xpd = xpd)
    on.exit(par(op))
    if (!add)
        plot(0, type = "n", ..., xlim = xlim, ylim = ylim, main = main,
            sub = sub, xlab = xlab, ylab = ylab, asp = 1, axes = axes)
    if (!plot)
        return()

    corm = cor(x)          #Distanzmatrix
    unm  = 1-abs(corm)
    unm[upper.tri(unm)]<-0
    hcl  = hclust(dist(unm))#clusterung
    bcl  = cutree(hcl,k=3)  #auf drei Cluster
    if(xA!=1 | yA!=2 | zA!=3) { 
        if(xA!=yA | yA!=zA | zA!=xA & (xA+yA+zA)==6) {
            bcl[bcl==1]=4   #Funktion zum tauschen
            bcl[bcl==2]=5   #auf gewünschte Reihenfolge
            bcl[bcl==3]=6
            bcl[bcl==4]=xA
            bcl[bcl==5]=yA
            bcl[bcl==6]=zA
        }
    }
    dmin=apply(x,2,min)     #Normierung des Datensatzes
    dmax=apply(x,2,max)
    q=(x-rep(1,nrow(x))%*%t(dmin))/(rep(1,nrow(x))%*%t(dmax-dmin))  
    points=matrix(0,nrow=3,ncol=2)
    for (i in 1:length(x[,1])) {
        xl=len*sum(q[i,bcl==1]) 
        yl=len*sum(q[i,bcl==2])
        zl=len*sum(q[i,bcl==3])
        points[1,1]=xloc[i]-yl/sqrt(2)
        points[1,2]=yloc[i]-yl/sqrt(2)
        points[2,1]=points[1,1]+xl
        points[2,2]=points[1,2]
        points[3,1]=xloc[i]+xl
        points[3,2]=yloc[i]
        lines(points,...)
    if (!is.null(labels)) { 
      text(points[1,1],min(points[,2]),labels[i],cex=cex,pos=1)
    }
        for (j in 1:length(q[1,bcl==3])) {
            c=q[i,bcl==3]
            points[,2]=points[,2]+len*as.numeric(c[j])
            lines(points,...) #Zeichnen der Merkmale 
        }                     #auf der z-Achse
        points[1,1]=xloc[i]-yl/sqrt(2)
        points[1,2]=yloc[i]-yl/sqrt(2)
        points[2,1]=points[1,1]
        points[2,2]=points[1,2]+zl
        points[3,1]=xloc[i]
        points[3,2]=yloc[i]+zl
        lines(points,...)
        for (j in 1:length(q[1,bcl==1])) {
            c=q[i,bcl==1]
            points[,1]=points[,1]+len*as.numeric(c[j])
            lines(points,...) #Zeichen der Merkmale
        }                     #auf der x-Achse
        points[1,1]=xloc[i]-yl/sqrt(2)+xl
        points[1,2]=yloc[i]-yl/sqrt(2)
        points[2,1]=points[1,1]
        points[2,2]=points[1,2]+zl
        points[3,1]=points[2,1]-xl
        points[3,2]=points[2,2]   
        for (j in 1:length(q[1,bcl==2])) {
            c=q[i,bcl==2]
            points[,1]=points[,1]+len*as.numeric(c[j])/sqrt(2)
            points[,2]=points[,2]+len*as.numeric(c[j])/sqrt(2)
            lines(points,...) #Zeichnen der Merkmale
        }                     #auf der y-Achse
    }
  

    if (!is.null(key.loc)) {
        par(xpd = key.xpd)

        q[1,]=1
        xleg=key.loc[1]
        yleg=key.loc[2]
        xl=leglen*sum(q[1,bcl==1])
        yl=leglen*sum(q[1,bcl==2])
        zl=leglen*sum(q[1,bcl==3])
        points[1,1]=xleg-yl/sqrt(2)
        points[1,2]=yleg-yl/sqrt(2)
        points[2,1]=points[1,1]+xl
        points[2,2]=points[1,2]
        points[3,1]=xleg+xl
        points[3,2]=yleg
        lines(points,...)
        for (j in 1:length(q[1,bcl==3])) {
            c=q[1,bcl==3]
            points[,2]=points[,2]+leglen*as.numeric(c[j])
            lines(points,...)
            text(xy.coords(points[1,1],points[1,2]),names(bcl[bcl==3])[j],cex=cex,
                 adj=c(1.3,1))
        }
        points[1,1]=xleg-yl/sqrt(2)
        points[1,2]=yleg-yl/sqrt(2)
        points[2,1]=points[1,1]
        points[2,2]=points[1,2]+zl
        points[3,1]=xleg
        points[3,2]=yleg+zl
        lines(points,...)
        for (j in 1:length(q[1,bcl==1])) {
            c=q[1,bcl==1]
            points[,1]=points[,1]+leglen*as.numeric(c[j])
            lines(points,...)
            text(xy.coords(points[1,1],points[1,2]),names(bcl[bcl==1])[j],cex=cex,
                 adj=c(1,1.5))
        }
        points[1,1]=xleg-yl/sqrt(2)+xl
        points[1,2]=yleg-yl/sqrt(2)
        points[2,1]=points[1,1]
        points[2,2]=points[1,2]+zl
        points[3,1]=points[2,1]-xl
        points[3,2]=points[2,2]   
        for (j in 1:length(q[1,bcl==2])) {
            c=q[1,bcl==2]
            points[,1]=points[,1]+leglen*as.numeric(c[j])/sqrt(2)
            points[,2]=points[,2]+leglen*as.numeric(c[j])/sqrt(2)
            lines(points,...)
            text(xy.coords(points[1,1],points[1,2]),names(bcl[bcl==2])[j],cex=cex,
                 adj=c(-0.2,1))
        }
    }
    if (frame.plot)
        box(...)
    invisible(locations)

}
