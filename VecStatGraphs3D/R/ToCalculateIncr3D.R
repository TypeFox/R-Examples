ToCalculateIncr3D <- function (coord) 
{
    num_data = dim(coord)
    num_elements = num_data[1]
    x_coordinates = 1:num_elements
    y_coordinates = 1:num_elements
    z_coordinates = 1:num_elements
    x_coordinates <- coord[, 4] - coord[, 1]
    y_coordinates <- coord[, 5] - coord[, 2]
    z_coordinates <- coord[, 6] - coord[, 3]
    increm = c(x_coordinates, y_coordinates, z_coordinates)
	dim(increm) = c(num_elements, 3)
    return(increm)
	
}
