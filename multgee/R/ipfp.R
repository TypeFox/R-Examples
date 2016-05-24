ipfp <-
function (initial, rowmars, colmars, dimension, maxit = ipfp.control()$maxit, 
    tol = ipfp.control()$tol) 
{
    ans <- initial
    rowsums <- .rowSums(ans, dimension, dimension, FALSE)
    for (i in 1:maxit) {
        ans <- ans * rep.int(rowmars/rowsums, dimension)
        colsums <- .colSums(ans, dimension, dimension, FALSE)
        if (all(abs(colsums - colmars) <= tol)) 
            break
        ans <- ans * rep(colmars/colsums, each = dimension)
        rowsums <- .rowSums(ans, dimension, dimension, FALSE)
        if (all(abs(rowsums - rowmars) <= tol)) 
            break
    }
    matrix(ans, dimension, dimension)
}

