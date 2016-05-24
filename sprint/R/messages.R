..sprintMsg <- list(
	error =
	c(non.double = "x must be of type double",
	non.square = "x is not and cannot be converted to a square matrix",
	non.ff = "x must be a valid ff object",
	non.ffmatrix = "ff object must be a matrix", 
	no.filename = "The filename of the ff object cannot be read",
	non.function = "The fun argument needs to be a function",
	non.supportedtype = "papply only supports a matrix or a list of matrices as input.",
	non.numeric = "PCOR only accepts numeric matrices",
	no.dims = "Dimensions of x and y matrices do not match",
	non.dna = "Function only accepts a character vector or an 'XStringSet' object",
	empty = "Data object is empty",
	no_file = "Output filename is missing",
	no.valid.k = "Number of clusters `k' must be in {1,2, .., n-1}; hence n >= 2",
	no.filename = "The filename of the ff object cannot be read",
	not.supported.non.hamming = "pstringdistmatrix only supports the hamming method. Please choose method=\"h\".",
	not.supported.diff.strings = "pstringdistmatrix only works when both input sets of strings are the same.",
	not.supported.maxDist = "maxDist is not used by pstringdistmatrix. Please set maxDist=0, or remove the maxDist parameter.",
	not.supported.ncores = "ncores is not used by pstringdistmatrix. Please refer to the SPRINT user guide for how to run in parallel."
	), 
	warn = c()
)