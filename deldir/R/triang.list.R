triang.list <- function (object) 
{
    stopifnot(inherits(object,"deldir"))
    io <- object$ind.orig
    pt <- object$summary$pt.type
    tlist <- triMat(object)
    x <- object$summary[,"x"]
    y <- object$summary[,"y"]
    if("z" %in% colnames(object$summary)) {
	z <- object$summary[,"z"]
	haveZ <- TRUE
    } else haveZ <- FALSE
    xtri <- matrix(x[tlist], nrow(tlist), 3)
    ytri <- matrix(y[tlist], nrow(tlist), 3)
    if(haveZ) ztri <- matrix(z[tlist], nrow(tlist), 3)
    ctri <- ytri - min(y)
    dx <- cbind(xtri[, 2] - xtri[, 1], xtri[, 3] - xtri[, 2], 
        xtri[, 1] - xtri[, 3])
    zm <- cbind(ctri[, 1] + ctri[, 2], ctri[, 2] + ctri[, 3], 
        ctri[, 3] + ctri[, 1])
    negareas <- apply(dx * zm, 1, sum)
    clockwise <- (negareas > 0)
    if (any(clockwise)) {
        xc <- xtri[clockwise,,drop=FALSE]
        yc <- ytri[clockwise,,drop=FALSE]
        tc <- tlist[clockwise,,drop=FALSE]
        if(haveZ) zc <- ztri[clockwise,,drop=FALSE]
        xtri[clockwise, ] <- xc[, c(1, 3, 2)]
        ytri[clockwise, ] <- yc[, c(1, 3, 2)]
        tlist[clockwise,] <- tc[, c(1, 3, 2)]
        if(haveZ) ztri[clockwise, ] <- zc[, c(1, 3, 2)]
    }
    rslt <- list()
    K <- 0
    for(i in 1:nrow(xtri)) {
	tmp <- .Fortran(
		"intri",
		x=as.double(xtri[i,]),
		y=as.double(ytri[i,]),
		u=as.double(x),
		v=as.double(y),
		n=as.integer(length(x)),
		okay=logical(1),
		PACKAGE="deldir"
	)
	if(tmp$okay) {
                tmp <- data.frame(ptNum=io[tlist[i,]],x=xtri[i,],y=ytri[i,])
                if(length(pt)) tmp <- cbind(tmp[,1,drop=FALSE],
                                            ptType=pt[tlist[i,]],tmp[,2:3])
                if(haveZ) tmp <- cbind(tmp,z=ztri[i,])
		K <- K+1
                rslt[[K]] <- tmp
        }
    }
    attr(rslt,"rw") <- object$rw
    class(rslt) <- "triang.list"
    rslt
}
"[.triang.list" <- function(x,i,...){
    y <- unclass(x)[i]
    class(y) <- "triang.list"
    attr(y,"rw") <- attr(x,"rw")
y
}
