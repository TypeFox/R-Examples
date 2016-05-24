hommel.test <-
#	Multiple testing from Hommel (1988).
#	Similar but very slightly more powerful that Hochberg (1988).
#	Controls Family-Wise Error rate for hypotheses which are independent or
#	which satisfy the free-association condition of Simes (1986).
#	Gordon Smyth, Walter and Eliza Hall Institute, smyth@wehi.edu.au
#	29 Aug 2002

function(p,alpha=0.05) {
	n <- length(p)
	i <- 1:n
	po <- sort(p)
	j <- n
	repeat {
		k <- 1:j
		if(all( po[n - j + k] > k * alpha / j )) break
		j <- j-1
		if(j == 0) break
	}
	p >= alpha/j
}
