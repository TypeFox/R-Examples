setmatch <-
function (x, y, nomatch = NA_integer_) 
{
    if (!is.list(x) && !is.list(y)) 
        stop("Arguments must be lists")
    out = rep.int(nomatch, length(x))
    for (i in seq_along(x)) {
        for (j in seq_along(y)) {
            if ((x[[i]] %subof% y[[j]]) && (y[[j]] %subof% x[[i]])) {
                out[i] = j
                break
            }
        }
    }
    out
}
