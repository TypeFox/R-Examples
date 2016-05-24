gromov.hyperbolicity <- function(d, deltas = FALSE, scale = NA)
{
    d <- as.matrix(as.dist(d))
    a = dim(d);
    if(a[1] != a[2])
    {
        stop("The parameter could not be coerced into a square distance matrix.")
    }

    if(a[1] < 4)
    {
        stop("At least 4 points must be used to compute the Gromov hyperbolicity constant.")
    }

    if(is.na(scale))
    {
        scale = "none"
    }

    scaleV = pmatch(scale, c("none", "max", "perimeter"))

    if(length(scaleV) > 1 || is.na(scaleV))
    {
        stop("Unknown or ambiguous scale method.");
    }

    .Call("gromov_distmatrix", as.double(d), as.logical(deltas),
            as.integer(scaleV), PACKAGE="distory")
}

