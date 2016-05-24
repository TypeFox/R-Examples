VectorsToPolar <- function (vectors) 
{
	
	num_data = dim(vectors)
    polar_vectors = vectors
    x = vectors[, 1]
    y = vectors[, 2]
    	
	module = sqrt(x * x + y * y)

    grades = ToSexagesimal(atan(x/y))
    gradesBool <- y >= 0

    grades[is.na(grades)] <- 0
    grades2 <- grades
    grades2[gradesBool == FALSE] <- grades[gradesBool == FALSE] + 
        180
    grades2Bool <- grades2 >= 0
    grades3 <- grades2
    grades3[grades2Bool == FALSE] <- grades2[grades2Bool == FALSE] + 
        360
    #polar_vectors = matrix(nrow = length(x), ncol = 2)
    polar_vectors[, 1] = module
    polar_vectors[, 2] = grades3
    return(polar_vectors)
}
