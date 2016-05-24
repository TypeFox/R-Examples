ToCalculateError <- function (matrix_) 
{
    num_data = dim(matrix_)
    num_elements = num_data[1]
    x_coordinates = 1:num_elements
    y_coordinates = 1:num_elements
    x_coordinates <- matrix_[, 1] - matrix_[, 3]
    y_coordinates <- matrix_[, 2] - matrix_[, 4]
    error = c(x_coordinates, y_coordinates)
    dim(error) = c(num_elements, 2)
    return(error)
}
