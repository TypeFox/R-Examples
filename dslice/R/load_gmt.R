load_gmt <- function(file = NULL)
{
	gmt_info <- readLines(file)
	set_num <- length(gmt_info)
	set_name <- vector(length = set_num, mode = "character")
	set_description <- vector(length = set_num, mode = "character")
	gene_symbol <- NULL
	for(i in 1:set_num){
		tmpset <- noquote(unlist(strsplit(gmt_info[i], "\t")))
		set_name[i] <- tmpset[1]
		set_description[i] <- tmpset[2]
		gene_symbol[[i]] <- tmpset[-(1:2)]
	}
	set_list <- list(set_name = set_name, set_description = set_description,
	                 gene_symbol = gene_symbol)
	return(set_list)
}
