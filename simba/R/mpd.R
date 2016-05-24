"mpd" <- 
function(x, method="simpson", all=FALSE)
{
	## preparations
	METHODS <- c("simpson", "sorensen", "nestedness")
	method <- pmatch(method, METHODS)
	## transform to presence/absence matrix
	df <- as.matrix(x > 0)
	## calculate b: species occurring only on one of any two compared plots
	b <- df %*% (1 - t(df))
	## calculate c: species occurring only on the other of any two compared plots
	c <- (1 - df) %*% t(df)
	## Summing up 'min(bi,bj)' for all possible pairings
	smin <- sum(pmin(b[row(b) > col(b)], c[row(c) > col(c)]))
	## Summing up 'max(bi,bj)' for all possible pairings
	smax <- sum(pmax(b[row(b) > col(b)], c[row(c) > col(c)]))
	## the number of species shared by at least two sites, with
	## sum(df) = sum of total species occurrences within data matrix,
	## 2nd term: total number of species including a check for empty species
	nss <- sum(df) - sum(colSums(df>0) != 0)
	
	## calculating multi plot dissimilarity indices
	sim <- smin / (nss + smin)
	sor <- (smin + smax) / (2*nss + smin + smax)
	nes <- sor-sim
	res <- c(sim, sor, nes)
	names(res) <- c("Simpson based mpd", "Sorensen based mpd", "Nestedness mpd")
	if(!all){
		res <- res[method]
	}
	res
}