#############################################################
#                                                           #
#   as.data.frame.circular function                         #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: September, 22, 2003                               #
#   Version: 0.1-2                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

as.data.frame.circular <- function(x, row.names=NULL, optional=FALSE, ...) {
    if (is.matrix(x)) {
        if (!is.null(xcircularp <- circularp(x))) {
            typep <- xcircularp$type
            unitsp <- xcircularp$units
            templatep <- xcircularp$template
            modulop <- xcircularp$modulo
            zerop <- xcircularp$zero
            rotationp <- xcircularp$rotation
        } else {
            typep <- "angles"
            unitsp <- "radians"
            templatep <- "none"
            modulop <- "asis"
            zerop <- 0
            rotationp <- "counter"
        }

        d <- dim(x)
        nrows <- d[1]; ir <- seq(length = nrows)
        ncols <- d[2]; ic <- seq(length = ncols)
        dn <- dimnames(x)
        row.names <- dn[[1]]
        collabs <- dn[[2]]
        if (any(empty <- nchar(collabs)==0))
        collabs[empty] <- paste("Circular", ic, sep = "")[empty]
        value <- vector("list", ncols)
    for(i in ic)
        value[[i]] <- as.circular(x[,i], type=typep, units=unitsp, modulo=modulop, zero=zerop, rotation=rotationp)
        if (length(row.names) != nrows)
        row.names <- if(optional) character(nrows) else as.character(ir)
        if (length(collabs) == ncols)
        names(value) <- collabs
        else if(!optional)
            names(value) <- paste("Circular", ic, sep="")
        attr(value, "row.names") <- row.names
        class(value) <- "data.frame"
        return(value)

    } else
    return(as.data.frame.vector(x, row.names, optional))
}
 
