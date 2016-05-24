subtable <-
function (x, variables, levels, drop = TRUE) 
{
    indexlist <- lapply(dim(x), seq_len)
    indexlist[variables] <- levels
    dims = dim(x)
    if (is.list(levels)) 
        dims[variables] = fsapply(levels, length)
    else dims[variables] = 1
    if (isTRUE(drop)) 
        dims = dims[!(seq_along(dims) %in% variables) | dims > 
            1]
    out = subarray(x, indexlist)
    dim(out) = dims
    return(out)
}
