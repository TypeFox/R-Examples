scanVectorMayNotExist <- function(fileName) {
	if (!file.exists(fileName)) {
		warning(paste("Cannot find file '", fileName, "'", sep=""))
		return(vector())
	}

	############################
	# AS IN XHMM C++ CODE:
	# Instead of splitting by whitespace, use ONLY tab as delimiter (to allow for sample names with space in them):
	############################

	return(scan(fileName, what=list(a=""), sep="\t")$a)
}
