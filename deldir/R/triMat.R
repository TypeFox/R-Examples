triMat <- function (object) {
#
# Function triMat to list the indices of the vertices of each
# Delaunay triangle in the triangulation of a planar point set.
# The indices are listed as the rows of an n x 3 matrix where n is
# the number of Delaunay triangles in the triangulation.
# 
    stopifnot(inherits(object, "deldir"))
    a  <- object$delsgs[, 5]
    b  <- object$delsgs[, 6]
    tlist <- matrix(integer(0), 0, 3)
    for (i in seq(nrow(object$summary))) {
        jj <- c(b[a == i], a[b == i])
        jj <- sort(unique(jj))
        jj <- jj[jj > i]
        if (length(jj) > 0)
            for (j in jj) {
                kk <- c(b[a == j], a[b == j])
                kk <- kk[(kk %in% jj) & (kk > j)]
                if (length(kk) > 0)
                  for (k in kk) tlist <- rbind(tlist, c(i, j, k))
            }
    }
    tlist
}
