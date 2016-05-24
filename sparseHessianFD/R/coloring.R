## Part of the sparseHessianFD package
## Copyright (C) 2013-2016 Michael Braun

#' @name coloring
#' @title Triangular partitioning of variables
#' @param L  sparsity pattern of the Hessian as a lower triangular pattern matrix
#' @return Integer vector of length nvars with color assignments for each variable.
#' @description cyclic coloring from a lower triangular pattern matrix
coloring <- function(L) {

    stopifnot(is(L,"nMatrix"),
              Matrix::isTriangular(L, upper=FALSE)
              )

    nvars <- NROW(L)
    L <- as(L,"nMatrix")
    G <- Matrix::crossprod(L)  # intersection graph
    ptr <- Matrix.to.Pointers(G, index1=FALSE)

    ## vertex coloring of intersection graph
    colors_vec <- get_colors(ptr[[2]], ptr[[1]], nvars)
    return(colors_vec)
}



