getfield <- function(data, field, whichline, pattern, replacement){
	tmp <- data[grep(field, data)+whichline] 
	tmp <- sub(tmp, pattern = pattern, replacement = replacement)
	return(tmp)
}