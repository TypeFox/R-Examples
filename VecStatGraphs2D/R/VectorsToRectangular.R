VectorsToRectangular <- function (vectors) 
{
    rectangular_vectors = vectors
    grades <- vectors[, 2]
    module <- vectors[, 1]
    radians <- ToRadians(grades)
    x1 = sin(radians) * module
    y1 = cos(radians) * module
	
	rectangular_vectors[, 1] <- x1
    rectangular_vectors[, 2] <- y1
    return(rectangular_vectors)
}
