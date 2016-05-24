cv.folds <-
function(n, folds = 10)
{
	split(sample(1:n), rep(1:folds, length = n))
}

