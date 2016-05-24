"image.khr" <-function(x, axes=FALSE, mar=c(0,0,2,0),
                       addcontour=TRUE, addpoints=TRUE,...)
{
    ## Verifications
    if (!inherits(x, "khr"))
        stop("x should be an object of class \"khr\"")
    if ((inherits(x,"khrud"))|(inherits(x,"kbbhrud")))
        col<-gray((256:1)/256)
    if (inherits(x,"khrvol"))
        col<-gray((1:256)/256)

    ## Graphical settings
    if (length(x) > 1) {
        opar<-par(mfrow=n2mfrow(length(x)), mar=mar)
        on.exit(par(opar))
    }

    ## For each animal, an image
    for (i in 1:length(x)) {
        if (length(x)>1)
            image(x[[i]]$UD, main=names(x)[i], axes=axes, col=col, ...)
        if (length(x)==1)
            image(x[[i]]$UD, axes=axes, col=col, ...)
        if (addcontour)
            contour(x[[i]]$UD, add=TRUE)
        if (addpoints) {
            if (inherits(x, "kbbhrud")) {
                points(x[[i]]$locs[[1]][,c("x","y")],
                       pch=21, col="black", bg="white")
            } else {
                points(x[[i]]$locs, pch=21, col="black", bg="white")
            }
        }
        box()
    }
}

