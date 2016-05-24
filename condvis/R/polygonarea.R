polygonarea <- function (x, y = NULL)
# calculate the area of polygon coming out of chull function
{
    if (is.null(y) && identical(ncol(x), 2L)){
        y <- x[, 2L]
        x <- x[, 1L]
    }
    area <- 0
    n <- length(x)
    j <- n
    for (i in 1:n){
        area <- area + (x[j] + x[i]) * (y[j] - y[i])
        j <- i
    }
    abs(area) / 2
}