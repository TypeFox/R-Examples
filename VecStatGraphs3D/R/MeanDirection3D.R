MeanDirection3D <- function (coord) 
{
	# the mean vector parameters are calculated by means 
	# of vectorial addition
	# unit vectors are not assumed
	
    x <- coord[, 1]
    y <- coord[, 2]
    z <- coord[, 3]
    n_elements = length(x)
    R = sqrt((sum(x) * sum(x)) + (sum(y) * sum(y)) + (sum(z) * 
        sum(z)))
    meanX <- sum(x) / R
    meanY <- sum(y) / R
    meanZ <- sum(z) / R
    if ((meanY > 0) & (meanX > 0)) {
        meanLongitud <- atan(meanY/meanX)
    }
    if ((meanY > 0) & (meanX < 0)) {
        meanLongitud <- atan(meanY/meanX) + pi
    }
    if ((meanY < 0) & (meanX < 0)) {
        meanLongitud <- atan(meanY/meanX) + pi
    }
    if ((meanY < 0) & (meanX > 0)) {
        meanLongitud <- atan(meanY/meanX) + 2 * pi
    }
    meanLongitud <- ToSexagesimal3D(meanLongitud)
    meanColatitud <- acos(meanZ)
    meanColatitud <- ToSexagesimal3D(meanColatitud)
    return(c(meanColatitud, meanLongitud))
}
