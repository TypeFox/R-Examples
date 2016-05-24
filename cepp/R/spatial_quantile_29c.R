#Core logic same as 29b, just rearranged to fit in better
evaluator <- function(n, p)
{
	phi <- function(u, t)	sqrt(sum(t^2)) + sum(u*t)
        
	I <- matrix(0, nrow=p, ncol=p)
	diag(I) <- 1
	I <- as.vector(I)

	func <- function(samp, at, u)
	{
		diff <- t(samp) - at
		scales <- sqrt(colSums(diff^2))
		scaled.diff <- t(diff) / scales

		result <- sum(scales + colSums(diff * u))
		
		delta <- colSums(scaled.diff) + n*u

		cross <- apply(scaled.diff, MARGIN=1, FUN=tcrossprod) #returns a column of values for each row of input
		inner <- I - cross
		sinner <- t(inner) / scales
		acc <- colSums(sinner)

		return(list(value=result, gradient=-delta, hessian=matrix(acc, ncol=p)))
	}
}
