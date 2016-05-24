"disjoin" <-
function(x) {
    # Transform factor levels into binary variables (one per level)
    if (!inherits(x, "factor"))
        stop("'x' must be a variable of class 'factor'!")
    N <- length(x)
	Nlevels <- nlevels(x)
	Xnum <- as.numeric(x)
	res <- matrix(rep(0, Nlevels * N), nrow = N)
	for (i in 1:Nlevels)
        res[Xnum == i, i] <- 1
    dimnames(res) <- list(names(x), levels(x))
    return(res)
}
