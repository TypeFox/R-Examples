
randomize.matrix <- function (x, ntimes=1, method=c("sample", "rowwise", "dataset", "complete"), seed=NULL, FUN=identity, ...) {
#-----------------------------------------------------------------------------------------
#  function leftover from an early idea for randomization tests,
#  currently unused but may become useful.  not exported.
#-----------------------------------------------------------------------------------------
	if(!is.matrix(x)) stop("argument must be matrix")
	f <- switch(
		match.arg (method),
		sample = function (x, ...) apply (x, 2, sample),
		rowwise = function (x, ...) apply (x, 1, sample),
		dataset = function (x, ...) matrix (sample (as.vector (x)), nrow (x), ncol (x)),
		complete = function (x, tot, ...)
				matrix (tabulate (sample (1:length(x), tot, TRUE), length (x)), nrow (x), ncol (x)))

	perms <- sapply (replicate (ntimes, x, simplify = FALSE), f, tot = sum (x), simplify = FALSE)
	sapply (perms, FUN, ..., simplify = FALSE)
	}
