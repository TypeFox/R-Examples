normalize <-
function(rawfeatures, meta=NULL) {
    if(is.null(meta) || is.null(meta$minValues)) {
        meta$minValues = apply(rawfeatures, 2, min, na.rm=T)
    }
    if(is.null(meta) || is.null(meta$maxValues)) {
        meta$maxValues = apply(rawfeatures, 2, max, na.rm=T)
    }
    tmp = rawfeatures
    for(i in names(rawfeatures)) {
        tmp[[i]] = if((meta$maxValues[[i]] - meta$minValues[[i]]) == 0) {
            rep.int(0, length(rawfeatures[[i]]))
        } else {
            ((rawfeatures[[i]] - meta$minValues[[i]]) / (meta$maxValues[[i]] - meta$minValues[[i]]) * 2) - 1
        }
    }

    return(list(features=tmp, meta=meta))
}
