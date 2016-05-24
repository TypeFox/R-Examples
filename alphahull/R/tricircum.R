tricircum <-
function (tri.obj) 
{
    if (!inherits(tri.obj, "tri")) 
        stop("tri.obj must be of class \"tri\"")
    nt <- summary(tri.obj)$nt
    ans <- .Fortran("delaunaycircum", as.integer(tri.obj$nc), 
        as.integer(tri.obj$lc), as.integer(tri.obj$n), as.double(tri.obj$x), 
        as.double(tri.obj$y), as.integer(tri.obj$tlist), as.integer(tri.obj$tlptr), 
        as.integer(tri.obj$tlend), as.integer(nt), circenter = double(2 * 
            nt), lct = integer(tri.obj$nc), tltri = integer(9 * 
            nt), ier = as.integer(0), PACKAGE = "alphahull")
    tri.info <- cbind(matrix(ans$tltri, nt, 9, byrow = TRUE), 
        matrix(ans$circenter, nt, 2, byrow = TRUE))
    colnames(tri.info) <- c("node1", "node2", "node3", "tr1", 
        "tr2", "tr3", "arc1", "arc2", "arc3", "circumx", "circumy")
    invisible(tri.info)
}
