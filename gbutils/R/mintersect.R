mintersect <- function(...){
    wrk <- list(...)
    if(length(wrk)==0) # 2015-02-12 was: return(character(0))
        stop("no arguments supplied, at least one is required")

    res <- wrk[[1]]
    if(length(wrk) > 1)
        for(i in 2:length(wrk))
            res <- intersect(res, wrk[[i]])
    res
}
