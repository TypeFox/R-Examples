downdateR <-
function(R, k = p)
{
	p <- dim(R)[1]
	if(p == 1)
		return(NULL)
	R <- delcol(R, rep(1, p), k)[[1]][ - p,  , drop = FALSE]
	attr(R, "rank") <- p - 1
	R	# Built-in Splus utility
}

