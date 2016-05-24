`as.array.mefa` <-
function(x, ...)
{
    DIM <- dim(x)
    DIMNAMES <- dimnames(x)
    DATA <- unlist(as.list.mefa(x))
    array(DATA, DIM, DIMNAMES)
}

