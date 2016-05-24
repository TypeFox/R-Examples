dichotomize.dataset <- function(x, split.at = NA) {
	if (is.na(split.at)) { split.at = median(x, na.rm = TRUE); }
	return( as.numeric( x > split.at ) );
	}


