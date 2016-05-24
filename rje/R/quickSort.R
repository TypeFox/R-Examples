quickSort <-
function (x, f = greaterThan, random = TRUE) 
{
    lst = is.list(x)
    n = length(x)
    if (n < 2) 
        return(seq_len(n))
    if (n == 2) {
        if (lst) 
            com = f(x[[1]], x[[2]])
        else com = f(x[1], x[2])
        if (com == 0) 
            return(c(1, 1))
        else if (com == 1) 
            return(c(2, 1))
        else if (com == -1) 
            return(c(1, 2))
        else stop()
    }
    if (random) 
        mid = sample(n, 1)
    else mid = ceiling(n/2)
    comp = numeric(n)
    comp[mid] = 2
    for (i in seq_len(n)[-mid]) {
        if (lst) 
            comp[i] = f(x[[i]], x[[mid]])
        else comp[i] = f(x[i], x[mid])
    }
    lu = Recall(x[comp == 1], f, random)
    ld = Recall(x[comp == -1], f, random)
    lm = Recall(x[comp == 0], f, random)
    rank = numeric(n)
    rank[comp == -1] = ld
    rank[mid] = max(c(0, ld)) + 1
    rank[comp == 0] = lm + max(c(0, ld))
    rank[comp == 1] = lu + max(c(0, rank), na.rm = TRUE)
    return(rank)
}
