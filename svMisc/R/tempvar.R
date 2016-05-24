tempvar <- function (pattern = ".var")
{
	## Similar to tempfile() but for temporary variables
	repeat {
		var <- paste(pattern, as.integer(runif(1) * 100000), sep = "")
		if (!exists(var, where = 1, inherits = TRUE)) break()
	}
	var
}
