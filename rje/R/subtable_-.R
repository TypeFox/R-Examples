`subtable<-` <-
function (x, variables, levels, value) 
{
    indexlist <- lapply(dim(x), seq_len)
    indexlist[variables] <- levels
    dims = dim(x)
    dims[variables] = sapply(levels, length)
    out = x
    subarray(out, indexlist) <- value
    return(out)
}
