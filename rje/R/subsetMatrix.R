subsetMatrix <-
function (n) 
{
    out = matrix(1, 1, 1)
    for (i in seq_len(n)) {
        out = matrix(c(1, -1, 0, 1), 2, 2) %x% out
    }
    out
}
