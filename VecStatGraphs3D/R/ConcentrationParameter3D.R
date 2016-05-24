ConcentrationParameter3D <- function (coord) 
{
    x <- coord[, 1]
    y <- coord[, 2]
    z <- coord[, 3]
    n_elements = length(x)
    mean_module = (MeanModule3D(coord))
    parameter = (n_elements - 1)/(n_elements * (1 - mean_module))
    return(parameter)
}
