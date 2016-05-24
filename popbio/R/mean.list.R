mean.list <- function(x, ...)
{
    if(!all(sapply(x, is.matrix)))
        stop("'x' must be a list containing matrices")
   dims <- sapply(x, dim)
   n <- dims[1, 1]
   p <- dims[2, 1]
   if(!all(n == dims[1, ]) || !all(p == dims[2, ]))
        stop("the matrices must have the same dimensions")
   mat <- matrix(unlist(x), n * p, length(x))
   mm<-matrix(rowMeans(mat, ...), n, p)
   dimnames(mm)<-dimnames(x[[1]])
   mm
}

