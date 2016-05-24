f.create.tag <- function(data, sep = "&"){
#
# CREATES A UNIQUE "TAG" FOR A COMBINATION OF COLUMNS BY
# PASTING INTO A SINGLE CHARACTER VECTOR
# NOTE: ORIGINAL DATA SHOULD not CONTAIN "&"S (THE SEPARATOR) FROM BEFORE...
#
# data CAN BE A data.frame OR A matrix
#
	.data.aux <- f.matrix.to.list(data)
	.data.aux <- lapply(.data.aux, as.character) # CONVERTS FACTORS ETC. TO CHARACTERS
	.data.aux$sep <- sep
	.tag <- do.call("paste", .data.aux)
	
	return(.tag)
}
