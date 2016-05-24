"extract.abund" <-
function(e, n, left=TRUE, ...) {
	if (missing(n)) n <- Inf
	if (is.null(n))
		stop("You must provide a value n for extraction!")
	p <- length(e$p.log.ind)
	# Verify n
	if (n < 0) n <- 0
	if (n > p) n <- p
	if (left == TRUE) { 			# Extract left part
		Res <- as.data.frame(eval(parse(text=e$data))[,e$vr][,1:n])
	} else {						# Extract right part
		Res <- as.data.frame(eval(parse(text=e$data))[,e$vr][,-(1:n)])
	}
	Res
}
