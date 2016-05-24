patternRepeat <-
function (x, which, n, careful = TRUE) 
{
    if (careful) {
        tmp = unique.default(which)
        if (!all(tmp == which)) {
            warning("repeated indices ignored")
            which = tmp
        }
        if (length(which) == 0) {
            if (length(x) != 1) 
                stop("x not of correct length")
            else return(rep.int(x, prod(n)))
        }
        if (max(which) > length(n)) 
            stop("n not long enough")
        if (prod(n[which]) != length(x)) 
            stop("x not of correct length")
    }
    else if (length(which) == 0) 
        return(rep.int(x, prod(n)))
    tmp = seq_along(n)[-which]
    if (length(tmp) == 0) 
        return(x)
    for (add.in in tmp) {
        bl = prod(n[seq_len(add.in - 1)])
        x = x[rep.int(seq_len(bl), n[add.in]) + rep(seq.int(from = 0, 
            by = bl, length.out = length(x)/bl), each = bl * 
            n[add.in])]
    }
    return(x)
}
