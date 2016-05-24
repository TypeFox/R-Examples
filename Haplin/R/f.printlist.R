f.printlist <- function(x){
##
##

.x <- lapply(x, function(y){
	if(is.list(y))return(y)
	else(return(list(y)))
})

.x <- substring(.x, first = 6, last = nchar(.x) - 1)
.x <- paste(names(x), ": ", .x, sep = "")

cat(.x, sep = "\n")

return(invisible(.x))

}
