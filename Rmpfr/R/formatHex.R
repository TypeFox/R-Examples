## sprintf("%+13.13a", x) ## hex digits after the hex point = 13

## precBits: double precision = 53 = 1 + 13*4


## conversion from Hex digits to binary sequences of digits
HextoBin <- c(
 "0"="0000",
 "1"="0001",
 "2"="0010",
 "3"="0011",
 "4"="0100",
 "5"="0101",
 "6"="0110",
 "7"="0111",
 "8"="1000",
 "9"="1001",
 "A"="1010",
 "B"="1011",
 "C"="1100",
 "D"="1101",
 "E"="1110",
 "F"="1111",
 "a"="1010",
 "b"="1011",
 "c"="1100",
 "d"="1101",
 "e"="1110",
 "f"="1111")

if(FALSE) {
## the code isn't using either of these inverses.
BintoHex <- names( HextoBin[1:16])
names(BintoHex) <- HextoBin[1:16]

Bintohex <- tolower(BintoHex)
}
##' @title sprintf("%a", *)-like formatting of mpfr numbers
##' @param x mpfr-number vector
##' @param bits integer (scalar) specifing the desired number of bits ("binary digits")
##' @param style 1-character string specifying
##' @return character vector of same length as \code{x}
##' @author Martin Maechler
sprintfMpfr <- function(x, bits, style = "+") {
    stopifnot(length(style <- as.character(style)) == 1, nchar(style) == 1,
	      style %in% c("+", " "),
	      length(bits) == 1, bits %% 1 == 0)
    if(bits > 52) {
	digits <- ceiling(bits / log2(10))
	neg <- sign(x) == -1
	ff <- .mpfr2str(x, digits, base = 16)
	isNum <- ff$finite	## ff$finite == is.finite(x)
	i0 <- ff$is.0	## == mpfrIs0(x)
	ex <- ff$exp ## the *decimal* exp : one too large *unless* x == 0
	r  <- ff$str # the mantissa, including "-" if negative
	Ex <- ex - 1L
	if(any(i0)) Ex[i0] <- ex[i0]
	if(!all(isNum)) ## "@Inf@", "@NaN@", ...
	    r[!isNum] <- gsub("@", '', r[!isNum], fixed=TRUE)
	r[isNum] <- paste0(substr(r[isNum], 1L, 1L), ".",
			   substring(r[isNum], 2L), "p", 4*Ex)
	if(any(i <- neg & isNum))
	    r[i] <- sub("^-", "-0x", r[i])
	if(any(i <- !neg & isNum))
	    r[i] <- paste0(style, "0x", r[i])
	r
    }
    else {
	hexdigits <- 1L + ((bits-1L) %/% 4L)
	nX <- as.character(hexdigits)
	sprintf(paste0("%", style, nX, ".", nX, "a"), x)
    }
}

formatHexInternal <- function(x, precBits = min(getPrec(x)), style = "+")
{
    bindigits <- as.integer(precBits) - 1L
    hexdigits <- 1L + ((bindigits-1L) %/% 4L)
    ## hexdigits is the number of hex digits after the precision point
    structure(sprintfMpfr(x, bits=bindigits, style=style),
              ##---------
	      bindigits = bindigits,
	      hexdigits = hexdigits)
}# formatHexInternal

##___ ../man/formatHex.Rd ___
##           ~~~~~~~~~~~~
formatHex <- function(x, precBits = min(getPrec(x)), style = "+") {
    structure(formatHexInternal(x, precBits=precBits, style=style),
	      dim = dim(x), dimnames = dimnames(x), class = "Hcharacter")
}

formatBin <- function(x, precBits=min(getPrec(x)), scientific = TRUE,
		      left.pad = "_", right.pad = left.pad, style = "+")
{
    H <- formatHexInternal(x, precBits=precBits, style=style)
    ## bindigits is number of binary digits after the precision point
    bindigits <- attr(H, "bindigits")
    hexdigits <- attr(H, "hexdigits")# *must* be correct = #{pure digits between "." and "p"}
    attributes(H) <- NULL
    S <- substr(H, 1, 1) # sign
    A <- substr(H, 4, 4)
    B <- substr(H, 6, 6+(hexdigits-1))
    ## assumes *always* an exponent "p" which is correct
    pow <- substr(H, 6+hexdigits+1, 1000000L)
    sB <- strsplit(B, "")
    rsB <- do.call(rbind, sB)
    hrsB <- HextoBin[rsB]
    dim(hrsB) <- dim(rsB)
    hrsBa <- apply(hrsB, 1, paste, collapse="")
    hrsBb <- substr(hrsBa, 1, bindigits)
    ## While this is a truncation,
    ## the mpfr conversion assures that
    ## only zero characters are truncated.
    if (!scientific) {
	powers <- as.integer(pow)
        Left <- -powers + max(powers, 2-precBits)
        Right <- powers - min(powers, precBits-1)
	D <- cbind(S, "0b", strrep(left.pad, Left),
		   A, hrsBb, strrep(right.pad, Right))
	D2 <- apply(D, 1, function(x) do.call(paste, list(x, collapse="")))
	ilft <- as.integer(max(Left) + min(powers)) + 4L
	res <- paste0(substr(D2,      1L,   ilft  ), ".",
		      substr(D2, ilft+1L, 1000000L))
    }
    else {
	res <- cbind(S, "0b", A, ".", hrsBb, "p", pow)
	res <- apply(res, 1, function(x) do.call(paste, list(x, collapse="")))
    }
    structure(res, dim = dim(x), dimnames = dimnames(x), class = "Bcharacter")
}

print.Bcharacter <- function(x, ...) {
    print(unclass(x), quote=FALSE, right=TRUE, ...)
    invisible(x)
}
## print.Hcharacter not needed: using S3method(*, *, *) trick in ../NAMESPACE !


## formatDec1 <- function(x, ...) {
##   H <- formatHexInternal(x, ...)
##   precBits <- attr(H, "bindigits") + 1
##   decdigits <- ceiling(log(2^precBits, 10)+1)
##   result <- format(x, digits=decdigits, ...)
##   class(result) <- class(H)
##   result
## }

## this is the version I sent to Martin
## formatDec2 <- function(x, ..., digits, nsmall) {
##   H <- formatHexInternal(x, ...)
##   precBits <- attr(H, "bindigits") + 1
##   decdigits <- ceiling(log(2^precBits, 10))
##   Hx <- scan(text=H, quiet=TRUE)
##   absHx <- abs(Hx)
##   minabsHx <- min(absHx[absHx != 0]) ## actual 0 is wanted here
##   logminabsHx <- log(minabsHx %% 1, 10)
##   nsmall <- if (is.infinite(logminabsHx))
##               (decdigits-1) - log(min(absHx[absHx != 0]), 10)
##             else
##               ceiling(-logminabsHx)
##   result <- format(Hx, digits=decdigits, nsmall=nsmall, ...)
##   dim(result) <- dim(x)
##   dimnames(result) <- dimnames(x)
##   class(result) <- class(H)
##   result
## }

### FIXME: also needs sprintf() --- *AND*  print(x,  digits= ...) is more flexible, works for all mpfr
### ------ This function is not needed, deprecate !
formatDec <- function(x, precBits = min(getPrec(x)), digits=decdigits,
                      nsmall=NULL, scientific=FALSE, style="+", ...) {
    H <- formatHexInternal(x, precBits=precBits, style=style)
    ## precBits <- attr(H, "bindigits") + 1
    decdigits <- ceiling(log(2^precBits, 10)) + 1
    Hx <- scan(text=H, quiet=TRUE) # <- FIXME ??? parse(.) or just  as.numeric() ?
    Hx.range <- range(abs(Hx[Hx != 0])) ## excluding hard zero
    if(!missing(digits))
	digits <- max(digits, decdigits)
    res <-
	if (scientific) {
	    sprintf(paste0("%", style, 1, ".", as.character(digits), "e"), Hx)
	}
	else {
	    if(is.null(nsmall)) nsmall <- -floor(log(Hx.range[1], 10))
	    if (nsmall <= 0) nsmall <- max(nsmall + digits, 0) ## effective only with scientific=FALSE
	    format(Hx, digits=digits, nsmall=nsmall, scientific=FALSE, ...)
	}
    structure(res,
	      dim = dim(x),
	      dimnames = dimnames(x),
	      class = "noquote")
}

## Attempted sprintf.  doesn't align with "f", so can't use here.
## formatDec3 <- function(x, ..., displaydigits=decdigits, fmttype=c("e","f")) {
##   H <- formatHexInternal(x, ...)
##   precBits <- attr(H, "bindigits") + 1
##   decdigits <- ceiling(log(2^precBits, 10))
##   Hx <- scan(text=H, quiet=TRUE)
##   absHx <- abs(Hx)

##   Hx.range <- range(abs(Hx[Hx != 0]))
##   fieldwidth <- ceiling(diff(log(Hx.range, 10)))
##   fmttype <- match.arg(fmttype)
##   precision <- log(Hx.range[1], 10)
##   precision <- if (precision >= 0) 0 else -floor(precision)
##   digits <- max(decdigits, displaydigits)
##   switch(fmttype,
##          "e"=format <- paste("%+", as.character(digits), ".", as.character(digits), "e", sep=""),
##          "f"=format <- paste("%+", as.character(digits), ".", as.character(precision), "f", sep="")
##          )

##   result <- sprintf(format, Hx)

##   dim(result) <- dim(x)
##   dimnames(result) <- dimnames(x)
##   class(result) <- class(H)
##   result
## }


mpfr.Bcharacter <- function(x, precBits, scientific = NA, ...) {
    ## was scanBin()
    if (is.na(scientific)) ## we look for a "p" exponent..
        scientific <- any(grepl("p", x, fixed=TRUE))
    class(x) <- NULL
    noPrec <- (missing(precBits) || is.null(precBits))
    if (!scientific) {
        if (noPrec) {
	    ## Need to count the number of actual bits, and do it before converting "_" to "0".
	    ## Exclude  Inf, NaN, NA, ..
	    x. <- x[x %in% c("NaN", "NA", "Inf", "+Inf", "-Inf")]
	    ## "Unpad" {my first(!) example where a pipe is slightly more elegant}
	    x. <- gsub("_", "",
		       sub("_+$", "",
			   sub("^[-+]?0b", "",
			       sub(".", "", x, fixed=TRUE))), fixed=TRUE)
	    precBits <- max(nchar(x.), 1)
        }
        x <- gsub("_", "0", x) ## TODO: chartr(.......)
    }
    else if (noPrec) { ## scientific -- find precBits from looking at string:
	## assume a format such as "+0b1.10010011101p+1"
	no.p <- -1L == (i.p <- as.vector(regexpr("p", x, fixed=TRUE)))
	if(any(no.p)) ## no "p" - infer precision from others
	    i.p[no.p] <- round(mean(i.p[!no.p]))
	if(length(unique(i.p)) == 1L) ## all are the same
	    i.p <- i.p[1]
	precBits <- i.p - 5L
    }
    mpfr(x, base = 2, precBits=precBits, ...)
}

## A mpfr() method for "Hcharacter" .. but not quite -- called from mpfr() when appropriate
mpfr.Hcharacter <- function(x, precBits, ...) {
    class(x) <- NULL
    if (missing(precBits) || is.null(precBits)) {
	## assume a format such as "+0x1.745d17ap-4"
	no.p <- -1L == (i.p <- as.vector(regexpr("p", x, fixed=TRUE)))
	if(any(no.p))
	    i.p[no.p] <- round(mean(i.p[!no.p]))
	if( all(duplicated(i.p)[-1])) ## all are the same
	    i.p <- i.p[1]
	precBits <- 1 + (i.p - 6) * 4
    }
    mpfr(x, base = 16, precBits=precBits, ...)
}
