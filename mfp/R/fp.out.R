fp.out <- function(..., pos)
{
    #
    # Version 1.0.1    28Aug01
    #
    x <- list(...)
    nx <- length(x)
    tabs <- 0
    cat("\n")
    for(i in 1:nx) {
        wordi <- x[[i]]
        if(pos[i] > tabs) {
            cat(rep("\t", (pos[i] - tabs)))
            tabs <- pos[i]
        }
        wordi <- wordi[!is.na(wordi)]
        # deal with powers 
        if(length(wordi) > 0) {
            cat(wordi)
            if(length(wordi) > 1) {
                nwordi <- 1 + sum(nchar(wordi))
            }
            else nwordi <- nchar(wordi)
            if(nwordi > 7) {
                tabs <- tabs + 1
            }
        }
        else cat("\t")
    }
    invisible()
}
