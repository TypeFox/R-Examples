## adds percent size to each factor levels
## author: Gilbert Ritschard

group.p <- function(group, weights=NULL){
    ##if(!is.factor(group)){
        group <- factor(group)
    ##}
    if(is.null(weights)) weights <- rep(1, length(group))
    ctb <- xtabs(weights ~ group)
    cprop <- prop.table(ctb)
    cprop.prct <- as.character(round(100*cprop,1))
    levels(group) <- paste(levels(group)," (",cprop.prct,"%)", sep="")
    return(group)
}
