MeanModule3D <- function (coord) 
{
    x <- coord[, 1]
    y <- coord[, 2]
    z <- coord[, 3]
	
    n_elements = length(x)
    R = sqrt((sum(x) * sum(x)) + (sum(y) * sum(y)) + (sum(z) * 
        sum(z)))
    mean_module = R/n_elements
    return(mean_module)
}
