Macdonald <-
function(x,m) {

x <- sqrt(2)*x

#
# Make some space for the part-answer
#
value <- matrix(0, nrow = length(x), ncol = m)

l.zero.handler <- function(x, m, j){

	x0.ix <- x==0

	answer <- rep(0, length(x))

	answer[!x0.ix] <- (m-j)*log(2*sqrt(m)*abs(x[!x0.ix]))

	if (m != j)
		answer[x0.ix] <- rep(log(0), sum(x0.ix))
	else
		answer[x0.ix] <- rep(0, sum(x0.ix))

	return(answer)

	}

for (j in 1:m) {

	value[,j] <- exp( l.zero.handler(x=x,m=m,j=j) + lfactorial(m+j-2)
		- lfactorial(m-j) - lfactorial(j-1)
		+ log(m)/2 - (2*m-1)*log(2) - lfactorial(m-1) - sqrt(m)*abs(x))
#	value[,j] <- exp( (m-j)*log(2*sqrt(m)*abs(x)) + lfactorial(m+j-2)
#		- lfactorial(m-j) - lfactorial(j-1)
#		+ log(m)/2 - (2*m-1)*log(2) - lfactorial(m-1) - sqrt(m)*abs(x))
	}	

answer <- rowSums(value)

return(answer*sqrt(2))
}
