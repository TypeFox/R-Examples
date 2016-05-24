## MM:  the "x <- a - a" (== integer(n)) and then
## ---      "x <- x + .C(...)"
## is probably for safety, normalization of unsigned int
## -- and for NA/NaN handling -- but that is done too complicated here.
##--- and: If I really want pack bits into integers,
##--- I cannot really deal with NA-bits in such a way!
##---> I'd need extra structure to store NA locations...
##
## OTOH: instead of CLASSES & COPY, we should use .Call() !!

bitFlip <- function(a,bitWidth=32)
{
    .Call("bitFlip", a, bitWidth, PACKAGE = "bitops")
}


bitAnd <- function(a, b)
{
    .Call("bitAnd", a, b, PACKAGE = "bitops")
}

bitOr <- function(a, b)
{
    .Call("bitOr", a, b, PACKAGE = "bitops")
}

bitXor <- function(a, b)
{
    .Call("bitXor", a, b, PACKAGE = "bitops")
}


bitShiftL <- function(a, b)
{
    .Call("bitShiftL", a, b, PACKAGE = "bitops")
}

bitShiftR <- function(a, b)
{
    .Call("bitShiftR", a, b, PACKAGE = "bitops")
}


cksum <- function(a)
{
    x <- nchar(as.character(a))*0
    x <- x + .C("cksum",
                length(a), as.character(a),
                val = as.numeric(x),
		NAOK=TRUE,
		DUP=TRUE,
                PACKAGE= "bitops")$val
    x[is.na(a)] <- NA
    x
}
