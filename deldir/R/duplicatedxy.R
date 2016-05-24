duplicatedxy <- function(x,y) {
    if(is.list(x)) {
        if(all(!is.na(match(c('x','y'),names(x))))) {
            return(duplicated(as.data.frame(x)))
        }
        stop("Argument \"x\" is a list but lacks x and/or y components.\n")
    }
    duplicated(data.frame(x=x,y=y))
}
