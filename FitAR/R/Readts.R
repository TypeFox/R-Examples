`Readts` <-
function(file = "", freq = 1, start = 1,  VerboseQ=TRUE)
{
	cat(title <- scan(file, what = "", sep = "\n", n = 1))
	commentcount <- 1
	while("#" == substring(scan(file, what = "", sep = "\n", n = 1, skip = 
		commentcount), first = 1, last = 1)) {
		commentcount <- commentcount + 1
	}
	if (VerboseQ) {
		cat("\n start = ")
		start <- as.numeric(eval(parse(text = readline())))
	}
	if (VerboseQ) {
		cat("\n frequency = ")
		freq <- as.integer(readline())
	}
	x <- scan(file = file, skip = commentcount)
	zts <- ts(x, start = start, frequency = freq)
	i <- nchar(title)
	while(substring(title, first = i, last = i) == " ") {
		i <- i - 1
	}
	title2 <- substring(title, first = 1, last = i)
	title2
	attr(zts, "title") <- title2
	zts
}

