codeword <-
function(word, m) {
result <- sort(word,index.return=T)
		permutation <- result$ix-1
		sorteddata <- result$x
		sortedlist <- (1:m)-1

		number <- 1
		for (j in 0:(m-1))
		{
			
			idx <- which(sortedlist==permutation[j+1])[1]
			sortedlist <- sortedlist[-idx]
			number <- number + factorial(m-j-1)*(idx-1)

			
		}	

return(number);
}