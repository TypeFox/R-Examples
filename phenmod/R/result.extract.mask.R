result.extract.mask <- function(mask.grid, values){
	# mask values: values in the area of germany should have valid value OR value -9999
	germany.mask <- which(is.na(mask.grid$alt)==FALSE)
		
	# values of doy.pim
	values.id.to.substitute <- germany.mask[which(is.na(values[germany.mask])==TRUE)]
	values[values.id.to.substitute] <- rep(-9999, length(values[values.id.to.substitute]))

	return(values)
}