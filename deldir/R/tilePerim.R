tilePerim <- function(object,inclbdry=TRUE) {
    if(!inherits(object,"tile.list"))
        stop("Argument \"object\" must be of class \"tile.list\".\n")
    perims <- sapply(object,tilePerim0,inclbdry=inclbdry)
    list(perimeters=perims,totalPerim=sum(perims),meanPerim=mean(perims))
}
