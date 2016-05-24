parseRd <- function (file)
{
	rdfile <- readLines(file)
	index <- cumsum(regexpr("^\\\\[[:alpha:]]", rdfile) > 0)
	
	chunks <- lapply(unique(index), function (i) rdfile[index == i])
	names <- sapply(chunks, function (x) gsub("^\\\\(.+)\\{.*$", "\\1", x[1]))
	
	cs <- cumsum(c(0, nchar(rdfile)))
	offset.start <- as.integer(cs[sapply(1:max(index),
		function (i) which(index == i)[1])])
	offset.end <- as.integer(cs[1 + sapply(1:max(index),
		function (i) tail(which(index == i), 1))])    
	
	return(list(offset.start = offset.start, offset.end = offset.end, 
		chunks = chunks, names = names))
}
