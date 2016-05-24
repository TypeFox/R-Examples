patterns <-
function(m, r, total) 
{
    y <- expand.grid(rep(list(0:(r-1)), m))
		y <- y[apply(y, 1, function(z) sum(z) %in% total),]
		colnames(y) <- paste("y", 1:m, sep = "")
		return(as.matrix(y))
}
