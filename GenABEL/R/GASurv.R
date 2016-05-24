"GASurv" <-
function(fuptime,status) {
	if (!is.numeric(fuptime))
		stop("follow-up time (fuptime) argument should be numeric!")
	tmp <- try(test.type(status,"binomial"))
	if (is(tmp,"try-error"))
		stop("Something is wrong with status argument")
	tmp <- try(test.type(fuptime,"gaussian"))
	if (is(tmp,"try-error"))
		stop("Something is wrong with fuptime argument")
	if (length(fuptime) != length(status))
		stop("fuptime and status arguments length differ")
	out <- matrix(c(fuptime,status),ncol=2)
	class(out) <- "GASurv"
	out
}
