##' Plot/estimate surface
##'
##' @export
##' @aliases ksmooth2 surface
##' @param x formula or data
##' @param data data.frame
##' @param h bandwidth
##' @param xlab X label
##' @param ylab Y label
##' @param zlab Z label
##' @param gridsize grid size of kernel smoother
##' @param ... Additional arguments to graphics routine (persp3d or persp)
##' @examples
##' ksmooth2(rmvn(1e4,sigma=diag(2)*.5+.5),c(-3.5,3.5),h=1,
##'         rgl=FALSE,theta=30)
##' 
##' if (interactive()) {
##'     ksmooth2(rmvn(1e4,sigma=diag(2)*.5+.5),c(-3.5,3.5),h=1)
##'     ksmooth2(function(x,y) x^2+y^2, c(-20,20))
##'     ksmooth2(function(x,y) x^2+y^2, xlim=c(-5,5), ylim=c(0,10))
##' 
##'     f <- function(x,y) 1-sqrt(x^2+y^2)
##'     surface(f,xlim=c(-1,1),alpha=0.9,aspect=c(1,1,0.75))
##'     surface(f,xlim=c(-1,1),clut=heat.colors(128))
##'     ##play3d(spin3d(axis=c(0,0,1), rpm=8), duration=5)
##' }
##' 
##' if (interactive()) {
##'     surface(function(x) dmvn(x,sigma=diag(2)),c(-3,3),lit=FALSE,smooth=FALSE,box=FALSE,alpha=0.8)
##'     surface(function(x) dmvn(x,sigma=diag(2)),c(-3,3),box=FALSE,specular="black")##' 
##' }
##' 
##' if (!inherits(try(find.package("fields"),silent=TRUE),"try-error")) {
##'     f <- function(x,y) 1-sqrt(x^2+y^2)
##'     ksmooth2(f,c(-1,1),rgl=FALSE,image=fields::image.plot)
##' }
ksmooth2 <- function(x,data,h=NULL,xlab=NULL,ylab=NULL,zlab="",gridsize=rep(51L,2),...) {
    if (is.function(x)) {
        args <- c(list(f=x,h=h,xlab=xlab,ylab=ylab,zlab=zlab),list(...))
        if (is.null(args$xlim) && !missing(data)) {
            if (is.list(data)) {
                args$xlim <- data[[1]]
                args$ylim <- data[[2]]
            } else args$xlim <- data
        }
        return(do.call(surface,args))
    }
    if (inherits(x,"formula")) {
        x <- model.frame(x,data)
    }
    if (length(gridsize)==1) gridsize <- rep(gridsize,2)
    if (is.null(h)) h <- sd(as.matrix(x))*nrow(x)^(-1/5)
    est <- KernSmooth::bkde2D(x, bandwidth=h, gridsize=gridsize)
    if (is.null(xlab)) xlab <- names(x)[1]
    if (is.null(ylab)) ylab <- names(x)[2]
    surface(est$fhat, x=est$x1, y=est$x2, est$fhat,
            xlab=xlab, ylab=ylab, zlab=zlab, ...)
    return(invisible(est))
}


##' @export
surface <- function(f,xlim=c(0,1),ylim=xlim,n=rep(100,2),col,clut="gold",clut.center,x,y,rgl=TRUE,expand=0.5,nlevels=10,col.contour="black",contour=TRUE,persp=TRUE,image="image",...) {
    if (missing(x)) {
        if (length(n)==1) n <- rep(n,2)
        x <- seq(xlim[1],xlim[2],length.out=n[1])
        y <- seq(ylim[1],ylim[2],length.out=n[2])
    }
    if (is.function(f)) {
        xy <- as.matrix(expand.grid(x,y))
        if (inherits(try(f(c(x[1],y[1])),silent=TRUE),"try-error")) {            
            f <- matrix(f(xy[,1],xy[,2]),nrow=length(x),ncol=length(y))
        } else {
            val <- f(xy)
            if (length(val)<NROW(xy)) {
                f <- matrix(apply(xy,1,f),nrow=length(x),ncol=length(y))
            } else {            
                f <- matrix(val,nrow=length(x),ncol=length(y))
            }
        }
    }
    zrg <- range(f)
    zlen <- diff(zrg)
    if (length(clut)==1) {
        ncolour <- 128
        clut <- switch(clut,
                       topo=topo.colors(ncolour),
                       red=colorRampPalette(c("yellow","red"),bias=1)(ncolour),
                       blue=colorRampPalette(c("white","blue"),bias=1)(ncolour),
                       gold=colorRampPalette(c("white","goldenrod1"),bias=1)(ncolour),
                       france=colorRampPalette(c("blue","white","red"))(ncolour),
                       rainbow=rainbow(n=128,start=0,end=1),
                       heat=heat.colors(ncolour),
                       heatrev=rev(heat.colors(ncolour)),                       
                       colorRampPalette(c("white","goldenrod1"),bias=1)(ncolour)
                       )
    }
    ncolour <- length(clut)
    if (!rgl) {
        if (missing(col)) {
            nc <- ncol(f); nr <- nrow(f)
            facet <- f[-1,-1]+f[-1,-nc]+f[-nr,-1]+f[-nr,-nc]
            facetcol <- cut(facet, ncolour)
            col <- clut[facetcol]
        }
        dots <- list(...)
        parargs <- list()
        if (persp) {
            parargs$mar <- c(0.2,0,0,0)
            if (contour || !is.null(image)) parargs$mfrow=c(2,1)
        }
        op <- do.call(par,parargs)
        if (persp) graphics::persp(x=x,y=x,z=f, col=col, expand=expand, ...)
        if (contour | !is.null(image)) {
            par(mar=c(3,3,0.5,3)) ##c(bottom, left, top, right)
            if (!is.null(image)) {
                do.call(image,list(x=x,y=y,z=f,col=clut,xlim=range(x),ylim=range(y),xlab="",ylab=""))
            }
            if (contour) {
                args <- c(list(x=x,y=y,z=f,nlevels=nlevels,col=col.contour,add=!is.null(image)),dots)
                args <- args[names(formals(contour.default))]
                do.call("contour",args) ##nlevels=20
            }
        }
        suppressWarnings(par(op))
    } else {
        if (missing(col)) {
            if (!missing(clut.center)) {
                zmax <- max(abs(zrg))
                zrg <- c(-zmax,zmax)
                zlen <- 2*zmax
            }
            col <- clut[round((length(clut)-1)*(f-zrg[1])/zlen)+1]
        }
        rgl::persp3d(x, y, f, col=col,...)
    }
}
