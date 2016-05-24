slidingWindow <-
function(DNAbin, width, interval = 1){
	annote <- function(x){
		width <- width - 1
		y <- align[ , x:(x+width)]
		attr(y, "window") <- c(x, x+width)
		y
	}
	if(interval == "codons") interval <- 3
	align <- as.matrix(DNAbin)
	len <- dim(align)[[2]]
	iter <- seq(1, len-width, by = interval)
	li <- lapply(iter, function(x) annote(x))
	li
}

