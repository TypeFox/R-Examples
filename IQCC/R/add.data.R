add.data <- function(datum2, estat, T2II, n, j, m = NULL)
{
    b <- T2.2(datum2, estat, n)
    if(is.null(m) == FALSE)
        j <- j + m + 1
    points(j, b[1], pch = 16)
    c <- c(j, j - 1)
    d <- c(b[1], T2II[1])
    lines(c, d)
    T2II <<- b
}