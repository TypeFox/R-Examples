na.stinterp <- function(object, ...) UseMethod("na.stinterp")

# interpolates object along along which defaults to time(object)
# along has to be numeric, is otherwise coerced
na.stinterp.default <- function(object, along = time(object), na.rm = TRUE, ...)
{
	along <- as.numeric(along)
	na.stinterp.0 <- function(y) {
		na <- is.na(y)
		if(all(!na)) return(y)
		y[na] <- stinterp(along[!na], as.numeric(y[!na]), along[na], ...)$y
		return(y)
	}
        object[] <- if (length(dim(object)) == 0) na.stinterp.0(object)
        	else apply(object, 2, na.stinterp.0)
        if (na.rm) na.omit(object) else object
}
