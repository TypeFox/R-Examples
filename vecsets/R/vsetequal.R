vsetequal <-
function (x, y, multiple=TRUE) {
    x <- as.vector(x)
    y <- as.vector(y)
	if(!multiple) {
		all(c(match(x, y, 0L) > 0L, match(y, x, 0L) > 0L))
	} else {
	# Can get away with this 'cause set theory doesn't "allow" floats
		length(x) == length(y) && identical(sort(x), sort(y))
	}
}
