plot.rPower <- function(x, ...) {
	if (!inherits(x, "rPower")) 
		stop("Use only with \"rPower\" objects")
  barplot(table(x[[2]]) / length(x[[2]]))
}

