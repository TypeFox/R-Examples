ordered_values <-
function(n) {

	sums <- function(m) {
		s <- 0
		for (i in 0:(m-1)) {
			s <- s + (1/(n - i))
			}
		s
		}

	ordered_values <- rep(NA, n)
	for (j in 1:n) {
		ordered_values[j] <- sums(j)
		}
	return(ordered_values)
	}
