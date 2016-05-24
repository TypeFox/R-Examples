
na.StructTS <- function(object, ...) UseMethod("na.StructTS")

na.StructTS.ts <- function(object, ..., na.rm = FALSE, maxgap = Inf)
{
    na.StructTS.0 <- function(y) {
        yf <- y
		isna <- is.na(y)
		yf[isna] <- rowSums(tsSmooth(StructTS(y))[,-2])[isna]
        .fill_short_gaps(y, yf, maxgap = maxgap)
    }
    object[] <- if (length(dim(object)) == 0) na.StructTS.0(object)
                else apply(object, 2, na.StructTS.0)
    if (na.rm) na.trim(object, is.na = "all") else object
}

na.StructTS.zoo <- function(object, ..., na.rm = FALSE, maxgap = Inf) {
	z <- na.StructTS(as.ts(object), ..., na.rm = FALSE, maxgap = maxgap)
	z <- as.zoo(z)
	time(z) <- time(object)
    if (na.rm) na.trim(z, is.na = "all") else z
}

