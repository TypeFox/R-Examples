streamorder <- function(dem,out) {

# data preparation

nrow <- dim(dem)[1]
ncol <- dim(dem)[2]

# calling the function

  result <- .C("streamorder",
		PACKAGE = "topmodel",
		as.double(dem),
		result = double(nrow*ncol),
		as.integer(nrow),
		as.integer(ncol),
		as.integer(out[1]),
		as.integer(out[2]))$result

# formatting of the results

return(matrix(result, nrow=nrow))

}
