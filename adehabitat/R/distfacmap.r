"distfacmap" <- function(x)
{
    ## Verifications
    if (!inherits(x, "asc"))
        stop("x should be of class \"asc\"")
    if (attr(x, "type")!="factor")
        stop("x should be of type \"factor\"")

    ## Bases for the function
    xyc <- getXYcoords(x)
    xc <- rep(xyc$x, times=length(xyc$y))
    yc <- rep(xyc$y, each=length(xyc$x))
    xyc<-data.frame(x=xc,y=yc)
    lev <- as.numeric(levels(factor(c(x))))
    li <- list()

    ## For each level of the map:
    for (i in lev) {

        ## keeps only the coordinates of the pixels
        ## corresponding to this level
        tmp <- x
        tmp[x!=i] <- NA
        tmp[x==i] <- 1
        ptsoui <- xyc[!is.na(c(tmp)),]

        ## these objects are passed to a call to the C function "distxyr",
        ## which computes the distance of each pixel to the nearest pixel
        ## for which the level is i
        toto <- .C("distxyr", as.double(t(as.matrix(xyc))),
                   as.double(t(as.matrix(ptsoui))),
                   as.integer(nrow(xyc)), as.integer(nrow(ptsoui)),
                   double(nrow(xyc)), PACKAGE="adehabitat")
        li[[i]] <- toto[[5]]
    }

    ## output as a kasc object
    names(li) <- levels(x)
    ka <- as.kasc(list(x1=x))
    li <- as.data.frame(li)
    li <- getkascattr(ka,li)
    li <- setmask(li, x)
    return(li)
}

