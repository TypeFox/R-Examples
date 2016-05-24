diff.square.sum <-
function(row1, rowdiff)
{
	# Helper function to compute the Brier score.
	y <- row1 - rowdiff
	sum(y*y)
}

