loadNamedVectorNoHeaderMayNotExist <- function(fileName) {
	if (!file.exists(fileName)) {
		warning(paste("Cannot find file '", fileName, "'", sep=""))
		return(vector())
	}

	############################
	# AS IN XHMM C++ CODE:
	# Instead of splitting by whitespace, use ONLY tab as delimiter (to allow for sample names with space in them):
	############################

	name_val = scan(fileName, what=list(name="", val=0.0), sep="\t")
	val = name_val$val
	names(val) = name_val$name

	return(val)
}
