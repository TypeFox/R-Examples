discretevar <-
function (data, var, n, p) 
{
    idx <- sort(data[, var], index.return = TRUE)
    x <- idx$x
    cx <- data[idx$ix, p]
    midpoint <- midpoints1(x)
    pts=rep(0,n)
    disc <- .C("discrete", as.double(x), as.integer(cx), as.double(midpoint), 
        as.integer(n), puntos = as.double(pts), npart = integer(1),PACKAGE="dprep")
    np <- disc$npart
    cat("The number of partitions for var", var, "  is :", np, 
        "\n")
    cat("The cut points are: ")
    points <- sort(disc$puntos[1:(disc$npart - 1)])
    print(points)
    y <- c(np, points)
}
