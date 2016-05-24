img <- function(x,idx,col=list(gray.colors(10,1,0.2)),
                ylab="Item",xlab="Subject",lab=TRUE,
                border=1,rowcol=FALSE,plotfun=NULL,
                axis1=TRUE,axis2=TRUE,yaxs="r",xaxs="r",cex.axis=0.4,...) {
    x0 <- seq(nrow(x))
    y0 <- seq(ncol(x))
    image(x=x0,y=y0,as.matrix(x),col=col[[1]],axes=FALSE,ylab=ylab,xlab=xlab,xaxs=xaxs,yaxs=yaxs,...)
    if (axis1) {
        axis(1,at=seq(nrow(x)),lwd=0.5,cex.axis=cex.axis,las=3)
        if (lab) suppressWarnings(title("",xlab=xlab,...))
    }
    if (axis2) {
        axis(2,at=seq(ncol(x)),lwd=0.5,cex.axis=cex.axis,las=1)
        if (lab) suppressWarnings(title("",ylab=ylab,...))
    }
    if (!is.null(plotfun)) {
        plotfun(...)
    }
    if (!missing(idx)) {
        if (rowcol) {
            for (i in seq_len(length(idx)))
                image(x=x0,y=idx[[i]],as.matrix(x[,idx[[i]]]),col=col[[i]],add=TRUE,xaxs=xaxs,yaxs=yaxs,...)
        } else
            for (i in seq_len(length(idx)))
                image(x=idx[[i]],y=y0,as.matrix(x[idx[[i]],]),col=col[[i]],add=TRUE,xaxs=xaxs,yaxs=yaxs,...)
    }
}



##' Visualize categorical by group variable
##'
##' @title Organize several image calls (for visualizing categorical data)
##' @param x data.frame or matrix
##' @param group group variable
##' @param ncol number of columns in layout
##' @param byrow organize by row if TRUE
##' @param colorbar Add color bar
##' @param colorbar.space Space around color bar
##' @param label.offset label offset
##' @param order order
##' @param colorbar.border Add border around color bar
##' @param main Main title
##' @param rowcol switch rows and columns
##' @param plotfun Alternative plot function (instead of 'image')
##' @param axis1 Axis 1
##' @param axis2 Axis 2
##' @param mar Margins
##' @param col Colours
##' @param ... Additional arguments to lower level graphics functions
##' @author Klaus Holst
##' @examples
##' X <- matrix(rbinom(400,3,0.5),20)
##' group <- rep(1:4,each=5)
##' images(X,colorbar=0,zlim=c(0,3))
##' images(X,group=group,zlim=c(0,3))
##' \dontrun{
##' images(X,group=group,col=list(RColorBrewer::brewer.pal(4,"Purples"),
##'                                RColorBrewer::brewer.pal(4,"Greys"),
##'                                RColorBrewer::brewer.pal(4,"YlGn"),
##'                                RColorBrewer::brewer.pal(4,"PuBuGn")),colorbar=2,zlim=c(0,3))
##' }
##' images(list(X,X,X,X),group=group,zlim=c(0,3))
##' images(list(X,X,X,X),ncol=1,group=group,zlim=c(0,3))
##' images(list(X,X),group,axis2=c(FALSE,FALSE),axis1=c(FALSE,FALSE),
##'       mar=list(c(0,0,0,0),c(0,0,0,0)),yaxs="i",xaxs="i",zlim=c(0,3))
##' @export
images <- function(x,group,ncol=2,byrow=TRUE,colorbar=1,colorbar.space=0.1,label.offset=0.02,
                 order=TRUE,colorbar.border=0,main,rowcol=FALSE,plotfun=NULL,
                   axis1,axis2,mar,
                   col=list(c("#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5"),
                       c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D"),
                       c("#EDF8E9", "#BAE4B3", "#74C476", "#238B45"),
                       c("#FEEDDE", "#FDBE85", "#FD8D3C", "#D94701")),
                   ...) {
    if (is.data.frame(x) || is.matrix(x)) x <- list(x)
    K <- length(x)
    lout <- matrix(seq(K),ncol=ncol,byrow=byrow)
    hei <- rep(1,nrow(lout))/nrow(lout)
    wid <- rep(1,ncol)/ncol
    if (colorbar==1) {
        wid <- c(rep(1,ncol)/ncol*(1-colorbar.space),colorbar.space)
        lout <- cbind(lout,K+1)
    }
    if (colorbar==2) {
        hei <- c(rep(1,nrow(lout))/nrow(lout)*(1-colorbar.space),colorbar.space)
        lout <- rbind(lout,K+1)
    }
    if (missing(group)) {
        group <- rep(1,nrow(x[[1]]))
    }
    if (missing(main)) main <- rep("",K)
    if (!is.list(col)) col <- list(col)
    group <- factor(group)
    idxs <- lapply(levels(group), function(x) which(group==x))
    layout(lout,widths=wid,heights=hei)
    ##if (missing(mar)) par(mar=c(4,4,3,0))
    if (missing(axis2)) axis2 <- c(TRUE,rep(FALSE,K-1))
    if (missing(axis1)) axis1 <- rep(TRUE,K)
    for (i in seq(length(x))) {
##        if (!missing(mar)) par(mar=mar[[i]])
        img(x[[i]],idxs,col,axis2=axis2[i],axis1=axis1[i],main=main[i],rowcol=rowcol,plotfun=plotfun[[i]],...)
##        if (missing(mar)) par(mar=c(4,2,3,2))
    }
    G <- nlevels(group)
    M <- length(col[[1]])
    if (colorbar==1) {
        par(mar=c(0,0,0,2))
        plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
        for (i in seq(G)) {
            lava::colorbar(col[[i]],values=seq(M)-1,direction="horizontal",
                           y.range=c(1-i/(G+1),1-i/(G+1)+label.offset),
                           border=colorbar.border,x.range=c(0,1),srt=0,cex=0.6)
            text(0.5,1-i/(G+1)-label.offset, levels(group)[i])
        }
    }
    if (colorbar==2) {
        par(mar=c(0,0,0,0))
        plot.new(); plot.window(xlim=c(0,1),ylim=c(0,1))
        for (i in seq(G)) {
            xr <- c(1-i/(G+1),1-i/(G+1)+.1)-.1/2
            lava::colorbar(col[[i]],values=seq(M)-1,direction="horizontal",
                           x.range=xr,
                           border=colorbar.border,y.range=c(0.3,0.5),srt=0,cex=0.6)
            text(mean(xr),.1, levels(group)[i])
        }
    }
}
