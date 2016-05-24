# save.cm validates and saves a Community Matrix. It takes:
# CM: a potential Community Matrix
# file: a valid filename with path
save.cm <- function(CM, file = stop("'file' must be specified")) {

	validate.cm(CM)
	save(CM, file = file, ascii = FALSE, compress = TRUE)
		
	# end Save.CM()
	}
