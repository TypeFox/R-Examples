track.plugin.lru <- function(objs, inmem, envname) {
    # A simple least-recently-used discard policy.
    # This function treats each tracked environment independently.
    # A more sophistocated version (to come) will look at all
    # tracked environments at once, and store a variable
    # ".trackingCacheMark" in each tracking env, which says which
    # vars to keep.
    cumsum.ordered <- function(x, order) return(replace(x, order, cumsum(x[order])))
    # just work with the objects that are in memory
    imobjs <- objs[inmem,]
    keep <- inmem
    # cache.size is in Mb
    cache.size <- getOption("track.cache.size")
    if (is.null(cache.size) || is.na(as.numeric(cache.size))) {
        cache.size <- NA
        if (.Platform$OS.type=="windows")
            cache.size <- memory.limit()/6
        if (is.na(cache.size))
            if (.Machine$sizeof.pointer > 4)
                cache.size <- 2048
            else
                cache.size <- 256
        options(track.cache.size=cache.size)
    }
    max.size <- (2^20) * cache.size
    # get the order, youngest first (prioritize ones with cache='yes' or 'fixedyes')
    by.age <- order(imobjs[,"cache"]=="yes" | imobjs[,"cache"]=="fixedyes", imobjs[,"accessed"], decreasing=TRUE)
    # which ones can we definitely keep?
    keep.by.age <- cumsum.ordered(imobjs[,"size"], by.age) <= max.size
    # now, order the potential deletion candidates by size to see if we can keep some
    # smaller ones one we get rid of the bigs ones
    by.size <- order(imobjs[!keep.by.age,"size"])
    keep.by.size <- cumsum.ordered(imobjs[!keep.by.age,"size"], by.size) <= (max.size - sum(imobjs[keep.by.age,"size"]))
    keep.by.age[which(!keep.by.age)] <- keep.by.size
    return(replace(inmem, which(inmem), keep.by.age))
}
