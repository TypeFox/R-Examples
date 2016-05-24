VectorsToPolar3D <- function(incr) 
{
    num_data = dim(incr)
    x = incr[, 1]
    y = incr[, 2]
    z = incr[, 3]
	
    module2D = sqrt(x^2 + y^2)
    colatitud = atan(module2D / z)
    colatitud[is.na(colatitud)] <- 0
    colatitudBool <- colatitud >= 0
    colatitud[colatitudBool == FALSE] <- colatitud[colatitudBool == 
        FALSE] + pi
    colatitud = ToSexagesimal3D(colatitud)
    longitud = ToSexagesimal3D(atan(y / x))
    longitud[is.na(longitud)] <- 0
    gradesBool <- x >= 0
    longitud[gradesBool == FALSE] <- longitud[gradesBool == FALSE] + 
        180
    gradesBool <- longitud >= 0
    longitud[gradesBool == FALSE] <- longitud[gradesBool == FALSE] + 
        360
    
	polar_vectors = matrix(nrow = num_data, ncol = 3)
	polar_vectors[, 1] = sqrt(x^2 + y^2 + z^2)
    polar_vectors[, 2] = colatitud
    polar_vectors[, 3] = longitud
    
	return(polar_vectors)
}
