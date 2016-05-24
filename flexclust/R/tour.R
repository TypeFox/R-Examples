#
#  Copyright (C) 2008 Friedrich Leisch
#  $Id: tour.R 3 2013-06-12 10:06:43Z leisch $
#

setGeneric("randomTour",
function(object, ...) standardGeneric("randomTour"))

setMethod("randomTour", signature(object="ANY"),
function(object, ...)
{
    object <- as(object, "matrix")
    randomTour(object, ...)
})

setMethod("randomTour", signature(object="matrix"),
function(object, ...)
{
    randomTourMatrix(object, ...)
})

setMethod("randomTour", signature(object="flexclust"),
function(object, data=NULL, col=NULL, ...)
{
    if(is.null(col)) col <- rep(flxColors(color="full"), length=object@k)
    
    if(is.null(data))
        data <- getData(object, error=TRUE)

    randomTourMatrix(data, col=col[predict(object, data)], ...)
})


###**********************************************************

randomTourMatrix <- function(x, directions=10,
                             steps=100, sec=4, sleep = sec/steps,
                       axiscol=2, axislab=colnames(x),
                       center=NULL, radius=1, minradius=0.01, asp=1,
                       ...)
{
    eins <- function(x) x/sqrt(sum(x^2))
    d <- ncol(x)
    D <- diag(d)
    axiscol <- rep(axiscol, d)

    ## arrows() gibt warning fuer falls pfeil zu kurs
    minradius <- max(minradius, .Machine$double.eps^0.25)

    if(is.null(center))
        docenter <- TRUE
    else{
        docenter <- FALSE
        center <- rep(center, length=2)
        Z <- matrix(center, byrow=TRUE, nrow=d, ncol=2)
    }

    ## zwei zufaellige richtungen und deren differenz
    A1 <- eins(rnorm(d))
    A2 <- eins(rnorm(d))
    A3 <- A2-A1

    for(m in 1:directions){
        for(n in (0:(steps-1))/steps){
            B <- as.matrix(eins(A1 + n*A3))
            ## Householder matrix
            B <- D - 2 * B %*% t(B)
            plot(x %*% B, xlab="", ylab="", asp=asp, axes=FALSE, ...)

            if(!is.null(axiscol)){
                if(docenter){
                    pu = par("usr")
                    center=c((pu[1]+pu[2])/2, (pu[3]+pu[4])/2)
                    Z <- matrix(center, byrow=TRUE, nrow=d, ncol=2)
                }
                C <- (D %*% B)[,1:2]
                CZ1 <- Z+radius*C
                CZ2 <- Z+1.1*radius*C
                ok <- rowSums(C^2)>minradius
                owarn <- options("warn")
                arrows(center[1], center[2], CZ1[ok,1], CZ1[ok,2],
                       col=axiscol[ok])
                text(CZ2[ok,1], CZ2[ok,2], labels=axislab[ok], col=axiscol[ok])
            }
            Sys.sleep(sleep)
        }
        A1 <- A2
        A2 <- eins(rnorm(d))
        A3 <- A2-A1
    }
}

###**********************************************************

randomTourRennes <- function(object, x=NULL, directions=10,
                             steps=100, sec=4, sleep = sec/steps,
                             axiscol=2, axislab=colnames(x),
                             center=NULL, radius=1, minradius=0.01, asp=1,
                             ...)
{
    if(is.null(x))
        x <- flexclust:::getData(object, error=TRUE)
    
    col <- clusters(object)
    
    eins <- function(x) x/sqrt(sum(x^2))
    d <- ncol(x)
    D <- diag(d)
    axiscol <- rep(axiscol, d)

    ## arrows() gibt warning fuer falls pfeil zu kurs
    minradius <- max(minradius, .Machine$double.eps^0.25)

    if(is.null(center))
        docenter <- TRUE
    else{
        docenter <- FALSE
        center <- rep(center, length=2)
        Z <- matrix(center, byrow=TRUE, nrow=d, ncol=2)
    }

    ## zwei zufaellige richtungen und deren differenz
    A1 <- eins(rnorm(d))
    A2 <- eins(rnorm(d))
    A3 <- A2-A1

    K <- object@k
    sim <- (object@clsim+t(object@clsim))/2
    sim[sim>0] <- 2

    for(m in 1:directions){
        for(n in (0:(steps-1))/steps){
            B <- as.matrix(eins(A1 + n*A3))
            ## Householder matrix
            B <- D - 2 * B %*% t(B)
            plot(x %*% B, xlab="", ylab="", asp=asp, axes=FALSE,
                 col=col, pch=19, cex=0.4, ...)

            points(object@centers %*% B, pch=1, cex=5, lwd=2)
            text(object@centers %*% B, labels=1:K, cex=1.5)

            for(k in 1:(K-1)){
                for(m in k:K)
                    if(sim[k,m]>0)
                        lines(object@centers[c(k,m),] %*% B, lwd=sim[k,m])
            }

            
            if(!is.null(axiscol)){
                if(docenter){
                    pu = par("usr")
                    center=c((pu[1]+pu[2])/2, (pu[3]+pu[4])/2)
                    Z <- matrix(center, byrow=TRUE, nrow=d, ncol=2)
                }
                C <- (D %*% B)[,1:2]
                CZ1 <- Z+radius*C
                CZ2 <- Z+1.1*radius*C
                ok <- rowSums(C^2)>minradius
                owarn <- options("warn")
                arrows(center[1], center[2], CZ1[ok,1], CZ1[ok,2],
                       col=axiscol[ok])
                text(CZ2[ok,1], CZ2[ok,2], labels=axislab[ok], col=axiscol[ok])
            }
            Sys.sleep(sleep)
        }
        A1 <- A2
        A2 <- eins(rnorm(d))
        A3 <- A2-A1
    }
}
