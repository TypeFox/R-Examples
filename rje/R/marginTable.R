marginTable <-
function (x, margin = NULL, order = TRUE) 
{
    k = length(dim(x))
    if (length(margin) > 0) 
        rmv = seq_len(k)[-margin]
    else return(sum(x))
    dx = dim(x)
    if (!is.double(x)) 
        x = as.double(x)
    out2 = .C("marginTable", x, as.integer(dx), as.integer(k), 
        as.integer(rmv), as.integer(length(rmv)), package = "rje")[[1]]
    if (order) {
        sm = sort.int(margin)
        kpm = dx[sm]
        out = out2[seq_len(prod(kpm))]
        dim(out) = kpm
        if (any(sm != margin)) 
            out = aperm.default(out, rank(margin))
        dimnames(out) = dimnames(x)[sm]
    }
    else {
        kpm = dx[margin]
        out = out2[seq_len(prod(kpm))]
        dim(out) = kpm
        dimnames(out) = dimnames(x)[margin]
    }
    out
}
