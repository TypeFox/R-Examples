caller <- function(start, index, n, bases)
{
	p <- vector('list', bases+1)
	i <- 1
        p[[1]] <- start
	cat("Initial Value", index(start), '\n')
	while (i < bases)
	{
		p[[i+1]] <- search_geodesic(current=p[[i]], index=index, n=n[i], max.tries=10)
		if(is.null(p[[i+1]])) return(p[1:i])
		i <- i + 1
	}
	return(p)
}
