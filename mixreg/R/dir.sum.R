dir.sum <- function(...)
{
# Function dir.sum.  To construct the direct sum of an arbitrary
# number of matrices.
	x <- list(...)
	x <- x[!sapply(x, is.null)]
	if(length(x) == 0)
		return(NULL)
	repeat {
		if(length(x) == 1) {
			if(is.list(x[[1]]))
				x <- x[[1]]
			else return(x[[1]])
		}
		else break
	}
	a <- as.matrix(x[[1]])
	b <- do.call("Recall", x[-1])
	ma <- nrow(a)
	na <- ncol(a)
	mb <- nrow(b)
	nb <- ncol(b)
	m <- ma + mb
	n <- na + nb
	rslt <- matrix(0, m, n)
	rslt[1:ma, 1:na] <- a
	rslt[(ma + 1):m, (na + 1):n] <- b
	rslt
}
