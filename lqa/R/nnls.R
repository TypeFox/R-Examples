nnls <- function (cvec, nnls.y, nnls.x)
{
    sum ((drop (nnls.y - nnls.x %*% cvec))^2)
}

