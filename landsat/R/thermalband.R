thermalband <-
function(x, band)
{
## convert thermal band from DN to temperature

	# coefs: gain, bias, K1, K2 from Chander et al. 2009
	if(band == 6) band.coefs <- c(0.055376, 1.18, 607.76, 1260.56)
	if(band == 61) band.coefs <- c(0.067087, -0.07, 666.09, 1282.71)
	if(band == 62) band.coefs <- c(0.037205, 3.16, 666.09, 1282.7)

	results <- x
	x <- as.vector(as.matrix(x))

	# at-sensor radiance
	x <- x * band.coefs[1] + band.coefs[2]

	x <- band.coefs[4] / log(band.coefs[3]/x + 1)

    # return the same structure as the input values
    if(class(results) == "SpatialGridDataFrame")
        results@data[,1] <- x
    else if(is.data.frame(results))
        results <- data.frame(matrix(x, nrow=nrow(results), ncol=ncol(results)))
    else if(is.matrix(results))
        results <- matrix(x, nrow=nrow(results), ncol=ncol(results))
    else # return a vector 
        results <- x
    
    results
}

