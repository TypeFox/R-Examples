getCol <- function (x,warn=FALSE) {
if(!inherits(x,"tile.list"))
    stop("Argument \"x\" must be of class \"tile.list\".\n")
ccc <- unlist(sapply(x,function(u){u[["z"]]}))
if(is.null(ccc)) return(NA)
ccc <- try(apply(col2rgb(ccc, TRUE), 2,
             function(x){do.call(rgb, as.list(x/255))}),silent=TRUE)
if(inherits(ccc,"try-error")){
    if(warn) warn(paste("Cannot interpret the z-components of",
                        "argument \"x\" as colours.\n"))
    return(NA)
}
ccc
}
