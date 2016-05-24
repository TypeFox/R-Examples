##' bin function
##'
##' A function to convert decimal to binary
##'
##' @param n a non-negative integer
##' @return the binary representation stored in a vector.
##' @export

bin <- function(n){
	m <- n
	ans <- NA
	while (m!=0){
		i <- 0
		while (2^i<=m){
			i <- i + 1
		}
		i <- max(i-1,0)
		ans[i+1] <- 1
		m <- m - 2^i
	}
	ans[ans!=1|is.na(ans)] <- 0
	return(rev(ans))
}
