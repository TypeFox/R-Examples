bal.table.ps <- function(x, digits = digits){
	lapply(x$desc, function(x){return(round(x$bal.tab$results, digits))})
}