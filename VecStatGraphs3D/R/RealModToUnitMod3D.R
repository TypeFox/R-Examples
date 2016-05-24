RealModToUnitMod3D <- function(incr)
{
# transforms vectors of arbitrary module to unit vectors

unit_incr = incr
n_elements = dim(incr)[1]
polar_values <- VectorsToPolar3D(incr)

module = matrix(1, n_elements, 1)
colatitude = polar_values[, 2]
longitude  = polar_values[, 3]
u_vector = matrix(c(module, colatitude, longitude), nrow = n_elements, ncol = 3)
unit_incr <- VectorsToRectangular3D(u_vector)
return(unit_incr)

}
	