getColors <-
function(n) {
    hues <- seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}

drawBox <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    UseMethod("drawBox")
}

drawBox.list <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
        lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    x <- x$pts
    NextMethod("drawBox")
}

drawBox.default <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
        lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    rect(x[1], x[2], x[3], x[4], col=bg, border=fg, lty=lty, lwd=lwd)
    ctr <- c(x[1] + (x[3]-x[1]) / 2, x[2] + (x[4]-x[2]) / 2)
    points(ctr[1], ctr[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)
}

drawBox2 <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    UseMethod("drawBox2")
}

drawBox2.list <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
        lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    x <- x$pts
    NextMethod("drawBox2")
}

drawBox2.default <-
function(x, fg=par("fg"), bg=NA, colCtr=NA,
        lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    polygon(x, col=bg, border=fg, lty=lty, lwd=lwd)
    ctr <- x[1, ] + 0.5 * (x[3, ] - x[1, ])
    points(ctr[1], ctr[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)
}

drawCircle <-
function(x, radius, nv=100, fg=par("fg"), bg=NA,
         colCtr=NA, lty=par("lty"), lwd=par("lwd"),
         pch=par("pch"), cex=par("cex")) {
    UseMethod("drawCircle")
}

drawCircle.list <-
function(x, radius, nv=100, fg=par("fg"), bg=NA,
         colCtr=NA, lty=par("lty"), lwd=par("lwd"),
         pch=par("pch"), cex=par("cex")) {
    radius <- x$rad
    x      <- x$ctr
    NextMethod("drawCircle", radius=radius)
}

drawCircle.default <-
function(x, radius, nv=100, fg=par("fg"), bg=NA,
         colCtr=NA, lty=par("lty"), lwd=par("lwd"),
         pch=par("pch"), cex=par("cex")) {
    if(!is.numeric(x))      { stop("x must be numeric") }
    if(!is.vector(x))       { stop("x must be a vector") }
    if(length(x) != 2L)     { stop("x must have length 2") }
    if(!is.numeric(radius)) { stop("radius must be numeric") }

    angles <- seq(0, 2*pi, length.out=nv)
    circ   <- cbind(x[1] + radius*cos(angles), x[2] + radius*sin(angles))

    polygon(circ[-1, ], border=NA, col=bg)
    lines(circ, col=fg, lwd=lwd, lty=lty)
    points(x[1], x[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)
}

drawEllipse <-
function(x, shape, radius, nv=100, axes=FALSE,
         fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    UseMethod("drawEllipse")
}

drawEllipse.list <-
function(x, shape, radius, nv=100, axes=FALSE,
        fg=par("fg"), bg=NA, colCtr=NA,
        lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    if(missing(shape))  { shape  <- x$cov }
    if(missing(radius)) { radius <- x$magFac }
    x <- x$ctr
    NextMethod("drawEllipse", shape=shape, radius=radius)
}

drawEllipse.default <-
function(x, shape, radius=1, nv=100, axes=FALSE,
         fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    if(!is.numeric(x))        { stop("x must be numeric") }
    if(!is.vector(x))         { stop("x must be a vector") }
    if(length(x) != 2L)       { stop("x must have length two") }
    if(!is.matrix(shape))     { stop("shape must be a matrix") }
    if(!is.numeric(shape))    { stop("shape must be numeric") }
    if(any(dim(shape) != 2L)) { stop("shape must be a (2 x 2)-matrix") }
    if(!isTRUE(all.equal(shape, t(shape)))) {
        stop("shape must be symmetric")
    }

    CF     <- chol(shape, pivot=TRUE)      # Cholesky-factor
    CFord  <- order(attr(CF, "pivot"))
    angles <- seq(0, 2*pi, length.out=nv)  # angles in radians
    ell    <- radius * cbind(cos(angles), sin(angles)) %*% CF[ , CFord]  # ellipse
    ellCtr <- sweep(ell, 2, x, "+")        # move ellipse to center

    ## draw center, ellipse
    points(x[1], x[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)  # center
    polygon(ellCtr, border=fg, col=bg, lwd=lwd, lty=lty)       # ellipse

    ## draw axes
    if(axes) {
        eig    <- eigen(shape)
        eigScl <- eig$vectors %*% diag(radius * sqrt(eig$values))

        # matrix with scaled ellipse axes
        xMat <- rbind(x[1] + eigScl[1, ], x[1] - eigScl[1, ])
        yMat <- rbind(x[2] + eigScl[2, ], x[2] - eigScl[2, ])

        matlines(xMat, yMat, col=fg, lwd=lwd, lty=lty)
    }
}

drawEllSector <-
function(x, shape=diag(2), radius=1, sect0=0, sect1=90, rot=0, nv=100,
         fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    if(!is.numeric(x))        { stop("x must be numeric") }
    if(!is.vector(x))         { stop("x must be a vector") }
    if(length(x) != 2L)       { stop("x must have length two") }
    if(!is.matrix(shape))     { stop("shape must be a matrix") }
    if(!is.numeric(shape))    { stop("shape must be numeric") }
    if(any(dim(shape) != 2L)) { stop("shape must be a (2 x 2)-matrix") }
    if(!isTRUE(all.equal(shape, t(shape)))) {
        stop("shape must be symmetric")
    }
    if(sect0 >= sect1) { stop("sect0 must be smaller than sect 1") }

    CF      <- chol(shape, pivot=TRUE)      # Cholesky-factor
    CFord   <- order(attr(CF, "pivot"))
    angles  <- seq(sect0*pi/180, sect1*pi/180, length.out=nv)  # angles in radians
    sect    <- radius * cbind(cos(angles), sin(angles)) %*% CF[ , CFord]  # sector
    sect    <- rbind(sect, c(0, 0))

    ang  <- -rot*pi/180
    G    <- cbind(c(cos(ang), sin(ang)), c(-sin(ang), cos(ang)))
    sect <- sect %*% G

    sectCtr <- sweep(sect, 2, x, "+")       # move sector to center

    ## draw center, ellipse
    points(x[1], x[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)  # center
    polygon(sectCtr, border=fg, col=bg, lwd=lwd, lty=lty)      # ellipse
}

drawTriSector <-
function(x, radius=1, sect0=0, sect1=90, rot=0,
         fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex")) {
    if(!is.numeric(x))  { stop("x must be numeric") }
    if(!is.vector(x))   { stop("x must be a vector") }
    if(length(x) != 2L) { stop("x must have length two") }
    if(sect0 >= sect1)  { stop("sect0 must be smaller than sect 1") }

    angles <- c(sect0*pi/180, sect1*pi/180)  # angles in radians
    sect   <- rbind(c(0, 0),
                    c(radius, radius*tan(angles[1])),
                    c(radius, radius*tan(angles[2])),
                    c(0, 0))

    ang  <- -rot*pi/180
    G    <- cbind(c(cos(ang), sin(ang)), c(-sin(ang), cos(ang)))
    sect <- sect %*% G

    sectCtr <- sweep(sect, 2, x, "+")       # move sector to center

    ## draw center, triangle
    points(x[1], x[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)  # center
    polygon(sectCtr, border=fg, col=bg, lwd=lwd, lty=lty)      # triangle
}

## draw oval shape from DSU targets
drawDSUOval <-
function(x, shape=diag(2), radius=1, angle, h=0, rot=0, nv=100,
         fg=par("fg"), bg=NA, colCtr=NA,
         lty=par("lty"), lwd=par("lwd"), pch=par("pch"), cex=par("cex"),
         plot=TRUE) {
    if(!is.numeric(x))        { stop("x must be numeric") }
    if(!is.vector(x))         { stop("x must be a vector") }
    if(length(x) != 2L)       { stop("x must have length two") }
    if(!is.matrix(shape))     { stop("shape must be a matrix") }
    if(!is.numeric(shape))    { stop("shape must be numeric") }
    if(any(dim(shape) != 2L)) { stop("shape must be a (2 x 2)-matrix") }
    if(!isTRUE(all.equal(shape, t(shape)))) {
        stop("shape must be symmetric")
    }
    if(angle >= 90) { stop("angle must be < 90 degree") }
    angle <- (pi/180)*angle     # convert angle to radians

    ## ellipse sector - right upper quadrant
    CF    <- chol(shape, pivot=TRUE)      # Cholesky-factor
    CFord <- order(attr(CF, "pivot"))

    ## correct opening angle for vertical offset h
    angCorr <- atan(h/radius)
    ellAng1 <- seq(angle-angCorr, (pi-angle+angCorr), length.out=nv)  # angles in radians
    ellAng2 <- pi + ellAng1
    ellSec1 <- radius * cbind(cos(ellAng1), sin(ellAng1)) %*% CF[ , CFord]  # sector
    ellSec2 <- radius * cbind(cos(ellAng2), sin(ellAng2)) %*% CF[ , CFord]  # sector

    ## move up to h
    ellSec1 <- sweep(ellSec1, 2, c(0,  h), "+")
    ellSec2 <- sweep(ellSec2, 2, c(0, -h), "+")

    ## triangle sector
    angTri  <- c(-angle, angle)
    triSec1 <- rbind(c( radius, radius*tan(angTri[1])),
                     c( radius, radius*tan(angTri[2])))
    triSec2 <- rbind(c(-radius, radius*tan(angTri[2])),
                     c(-radius, radius*tan(angTri[1])))

    ## join triangle - ellipse - triangle - ellipse
    oval <- rbind(triSec1,
                  ellSec1,
                  triSec2,
                  ellSec2)

    ## rotate and move to center
    ang     <- -rot*pi/180
    G       <- cbind(c(cos(ang), sin(ang)), c(-sin(ang), cos(ang)))
    oval    <- oval %*% G
    ovalCtr <- sweep(oval, 2, x, "+")       # move sector to center

    ## draw center, triangle
    if(plot) {
        points(x[1], x[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)  # center
        polygon(ovalCtr, border=fg, col=bg, lwd=lwd, lty=lty)      # triangle
    } else {
        return(ovalCtr)
    }
}

## draw a super ellipse
drawSuperEll <-
function(x, a, b, n, nv=100, rot=0,
         col=par("fg"), colCtr=NA, lty=par("lty"), lwd=par("lwd"),
         pch=par("pch"), cex=par("cex")) {
    angles <- seq(0, pi/2, length.out=nv)
    x0 <- a*cos(angles)^(2/n)
    y0 <- b*sin(angles)^(2/n)
    xy <- cbind(c(x0, rev(-x0), -x0, rev( x0)),
                c(y0, rev( y0), -y0, rev(-y0)))

    ang <- -rot*pi/180
    G   <- cbind(c(cos(ang), sin(ang)), c(-sin(ang), cos(ang)))
    xy  <- xy %*% G

    xyCtr <- sweep(xy, 2, x, "+")       # move sector to center

    lines(xyCtr, col=col, lwd=lwd, lty=lty)
    ## draw center
    points(x[1], x[2], col=colCtr, pch=pch, lwd=lwd, cex=cex)  # center
}
