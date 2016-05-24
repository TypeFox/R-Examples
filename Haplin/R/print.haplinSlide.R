print.haplinSlide <- function(x, ...){

.names <- names(x)
.na <- is.na(x)

for(i in seq(along = x)){
	cat("-----------------------\n")
	cat("Result from Window '", .names[i], "'...:\n", sep = "")
	if(.na[i]) cat("RUN FAILED\n")
	else print(x[[i]])
}
cat("-----------------------\n")


if(any(.na)){
	.mess <- sapply(x[.na], function(x) attr(x, "error.message"))
	cat("\nNOTE: The following window runs failed:\n")
	for(i in seq(along = .mess)){
		cat("=========\n")
		cat("Window: '", names(.mess)[i], "', Error message:\n", sep = "")
		cat(.mess[[i]])
	}
	cat("=========\n")
	###.res.list[.tmp1] <- NULL
}




}
