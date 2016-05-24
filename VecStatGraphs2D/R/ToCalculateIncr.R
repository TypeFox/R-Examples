ToCalculateIncr <- function (matrix_) 
{
    num_data = dim(matrix_)
    num_elements = num_data[1]
    x_coordinates = 1:num_elements
    y_coordinates = 1:num_elements
    x_coordinates <- matrix_[, 3] - matrix_[, 1]
    y_coordinates <- matrix_[, 4] - matrix_[, 2]
    incr = c(x_coordinates, y_coordinates)
    dim(incr) = c(num_elements, 2)
    return(incr)
}
