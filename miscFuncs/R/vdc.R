##' vdc function
##'
##' A function to generate a Van der Corput sequence of numbers.
##'
##' @param n the length of the sequence
##' @return Van der Corput sequence of length n
##' @export

vdc <- function(n){
	if(n==1){
		return(1/2)
	}
	else{
		denom <- 2
		i <- 2
		count <- 1
		while(length(denom)<n){
			denom <- c(denom,2*i)
			count <- count + 1
			if (count > i){
				count <- 1
				i <- 2*i
			}
		}
		num <- denom / 2
		const <- 0
		N <- floor(log(n,base=2))
		for (i in 1:N){
			s <- seq(1,(2^i)-1,by=2)
			le <- length(s)
			if (le>1){
				s <- c(s[le],s[1:(le-1)])
			}
			s <- c(s,rev(s))
			idx <- 1:length(s)
			s[idx%%2==1] <- -s[idx%%2==1]
			const <- c(const,s) 
		}
		const <- const[1:n]
		num <- num + const
	}
	return(num/denom)
}
