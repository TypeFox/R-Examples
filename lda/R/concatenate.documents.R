concatenate.documents <-
function (...) 
{
    lengths <- sapply(list(...), length)
    stopifnot(all(lengths == lengths[1]))
    mapply(cbind, ...)
}
