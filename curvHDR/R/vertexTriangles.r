
vertexTriangles <- function (ve)
{
    n.vert <- ncol(ve$vb)
    val <- vector("list", n.vert)
    ib <- ve$ib
    for (i in 1:ncol(ib)) {
        val[[ib[1, i]]] <- c(val[[ib[1, i]]], i)
        val[[ib[2, i]]] <- c(val[[ib[2, i]]], i)
        val[[ib[3, i]]] <- c(val[[ib[3, i]]], i)
    }
    return(val)
}

