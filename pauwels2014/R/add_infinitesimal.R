add_infinitesimal <-
function(theta, coords, coords_minus = c(), step = 1e-6, ...){
	## Adds or deletes the small quantitiy step to vector theta
	## at some coordinates
	temp <- theta
	ucoord_plus <- unique(coords)
	mult_coord_plus <- sapply(ucoord_plus, FUN = function(x){sum(coords == x)})
	ucoord_minus <- unique(coords_minus)
	mult_coord_minus <- sapply(ucoord_minus, FUN = function(x){sum(coords_minus == x)})
	if(length(mult_coord_plus) == 0){
		mult_coord_plus <- NULL
	} 
	if(length(mult_coord_minus) == 0){
		mult_coord_minus <- NULL
	}
	
	temp[ucoord_plus] <- temp[ucoord_plus] + step * mult_coord_plus
	temp[ucoord_minus] <- temp[ucoord_minus] - step * mult_coord_minus
	temp
}
