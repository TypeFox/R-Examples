
fromto <- function (x, from, to)
{
    if (is.list(x)) {
	xattr <- attributes(x)
	x <- lapply(x, fromto, from, to)
	attributes(x) <- xattr
    } else {
	if (is.factor(x)) {
	    levels(x) <- iconv(levels(x), from, to, sub = "byte")
	} else {
	    if (is.character(x))
		x <- iconv(x, from, to, sub = "byte")
	}
	lb <- attr(x, "label")
	if (length(lb) > 0) {
	    attr(x, "label") <- iconv(attr(x, "label"), from, to, sub = "byte")
	}
    }
    x
}

# Converts a variable from UTF-8 into other encoding
fromUTF8 <- function (x, to = "WINDOWS-1252")
{
    fromto(x, "UTF-8", to)
}

# Converts a variable from any encoding into UTF-8
toUTF8 <- function (x, from = "WINDOWS-1252")
{
    fromto(x, from, "UTF-8")
}

