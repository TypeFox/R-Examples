"bcoov" <-
	function(x, names, listout=FALSE) {
	anz <- length(x)
	bc <- abs(outer(x, x, "-")) / outer(x, x, "+")
	d <- as.dist(bc)
	attr(d, "Size") <- anz
    attr(d, "Labels") <- names
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    if (listout) {
    	d <- liste(d, entry="Bray-Curtis")
    	}
    return(d)
}