new.seed <- function() {
	x <- format(Sys.time(), "%H %M %S")
	x <- as.numeric(substring(x, 1, 2)) * 24 * 60 + as.numeric(substring(x, 4, 5)) * 60 + as.numeric(substring(x, 7, 8))
	x
}