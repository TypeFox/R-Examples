#### All  coercion methods for the  "Rmpfr" classes

if(getRversion() < "2.15")
    paste0 <- function(...) paste(..., sep = '')
if(getRversion() < "3.3")
    strrep <- function (x, times) { ## (x, times) must be "recycled"
	if((lx <- length(x)) < (lt <- length(times)))
	    x <- rep_len(x, lt)
	else if(lt < lx)
	    times <- rep_len(times, lx)
	vapply(seq_along(x),
	       function(i) paste(rep.int(x[i], times[i]), collapse = ""), "")
    }

mpfr <- function(x, precBits, ...) UseMethod("mpfr")

mpfr.mpfr <- function(x, precBits, rnd.mode = c('N','D','U','Z','A'), ...)
    roundMpfr(x, precBits=precBits, rnd.mode=rnd.mode)

mpfr.bigz <- function(x, precBits, ...) {
    if(missing(precBits)) precBits <- max(2L, frexpZ(x)$exp)
    if(getOption("verbose"))
	warning("mpfr(<bigz>) --> .bigz2mpfr() [not efficiently via character]")
    ..bigz2mpfr(x, precBits)
}

mpfr.bigq <- function(x, precBits, ...) {
    if(missing(precBits)) precBits <- getPrec(x)#-> warning
    if(getOption("verbose"))
	warning("mpfr(<bigq>) --> .bigq2mpfr() [not efficiently via character]")
    ..bigq2mpfr(x, precBits)
}

mpfr.default <- function(x, precBits, base = 10, rnd.mode = c('N','D','U','Z','A'),
                         scientific = NA, ...)
{
    if(is.ch <- is.character(x))
	stopifnot(length(base) == 1, 2 <= base, base <= 62)
    stopifnot(is.character(rnd.mode <- toupper(rnd.mode)))
    rnd.mode <- match.arg(rnd.mode)

    if(is.raw(x)) { # is.raw() is faster
	stopifnot(missing(precBits) || precBits >= 2)
	## else warning("unrecognized raw 'x'") # <- ?? {see use in ../tests/create.R }
        ## {but 'raw' is treated below}
    } ## else

    if(missing(precBits)) {
	precBits <- getPrec(x, base = base, doNumeric = FALSE)
    }
    ## libmpfr would exit (after good error message) for precBits == 1
    stopifnot(precBits >= 2)

    ml <-
	if(is.numeric(x) || is.logical(x) || is.raw(x))
	    .Call(d2mpfr1_list, x, precBits, rnd.mode)
	else if(is.ch)
	    .Call(str2mpfr1_list,x, precBits, base, rnd.mode)
	else stop("invalid 'x'. Must be numeric (logical, raw) or character")
    if(is.array(x)) {
	dim <- dim(x) ; dn <- dimnames(x)
	new(if(length(dim) == 2) "mpfrMatrix" else "mpfrArray",
	    ml, Dim = dim,
	    Dimnames = if(is.null(dn)) vector("list", length(dim)) else dn)
    }
    else new("mpfr", ml)
} ## mpfr.default()

.mpfr <- function(x, precBits)
    new("mpfr", .Call(d2mpfr1_list, x, precBits, "N"))
.mpfr. <- function(x, precBits, rnd.mode)
    new("mpfr", .Call(d2mpfr1_list, x, precBits, rnd.mode))


##' to be used in our own low-level R programming
.d2mpfr1 <- function(x, precBits) .Call(d2mpfr1, x, precBits, "N")
setAs("numeric", "mpfr1", ## use default precision of 128 bits
      function(from) .Call(d2mpfr1, from, 128L, "N"))# <- round to [N]earest
setAs("numeric", "mpfr", function(from) .mpfr(from, 128L))
setAs("integer", "mpfr", function(from) .mpfr(from,  32L))
setAs("raw",     "mpfr", function(from) .mpfr(from,   8L))
setAs("logical", "mpfr", function(from) .mpfr(from,   2L))
## TODO?  base=16 for "0x" or "0X" prefix -- but base must have length 1 ..
setAs("character", "mpfr", function(from) mpfr(from))

setAs("mpfr", "numeric", function(from) .Call(mpfr2d, from, rnd.mode="N"))
setAs("mpfr", "integer", function(from) .Call(mpfr2i, from, rnd.mode="N"))
setMethod("as.numeric", "mpfr", function(x, rnd.mode="N") .Call(mpfr2d, x, rnd.mode))
## "Z": round towards [Z]ero -- crucial for as.integer() :
setMethod("as.integer", "mpfr", function(x, rnd.mode="Z") .Call(mpfr2i, x, rnd.mode))

## FIXME (in gmp!!): asNumeric() should get "..." argument
setMethod("asNumeric", "mpfr",      function(x) .Call(mpfr2d, x, rnd.mode="N"))
setMethod("asNumeric", "mpfrArray", function(x) toNum(x, rnd.mode="N"))

setAs("mpfr1", "numeric",  ## just for user-de-confusion :
      function(from) {
	  warning("coercing \"mpfr1\" via \"mpfr\" (inefficient)")
	  as(new("mpfr", list(from)), "numeric") })

setAs("mpfr1", "mpfr", function(from) new("mpfr", list(from)))
setAs("mpfr", "mpfr1", function(from) {
    if(length(from) == 1) getD(from)[[1]] else
    stop("only \"mpfr\" objects of length 1 can be coerced to \"mpfr1\"")
})

.mpfr1tolist <- function(x)
    sapply(.slotNames(x), slot, object=x, simplify=FALSE)
.mpfr2list <- function(x) lapply(getD(x), .mpfr1tolist)

## Breaks the working of vapply(q, FUN.x) in pbetaI() in ./special-fun.R :
## as.list.mpfr1 <- function(x, ...) .mpfr1tolist(x)
## as.list.mpfr  <- function(x, ...) .mpfr2list(x)

## and then
mpfrXport <- function(x) {
    if(!is(x, "mpfr")) stop("argument is not a \"mpfr\" object")
    structure(class = "mpfrXport",
	      list(gmp.numb.bits = .mpfr.gmp.numbbits(),
		   ## currently unused, but in case:
		   mpfr.version	 = .mpfrVersion(),
		   Machine  = .Machine[grepl("sizeof",names(.Machine))],
		   Sys.info = Sys.info()[c("sysname", "machine")],
		   mpfr = .mpfr2list(x)))
}

mpfrImport <- function(mxp) {
    if(!inherits(mxp, "mpfrXport")) stop("need an \"mpfrXport\" object")
    nbits <- .mpfr.gmp.numbbits()
    if(!identical(nbits, mxp$gmp.numb.bits))
	stop("GMP bits not matching: 'x' has ", mxp$gmp.numb.bits,
	     "; the loaded 'Rmpfr' package has ", nbits)
    m1 <- lapply(mxp$mpfr, function(o) do.call(new, c("mpfr1", o)))
    new("mpfr", m1)
}

.mpfr2str <- function(x, digits = NULL, base = 10L) {
    ##	digits = NULL : use as many digits "as needed"
    stopifnot(is.null(digits) ||
	      (is.numeric(digits) && digits >= 1),
	      is.numeric(base), length(base) == 1, base >= 2)
    .Call(mpfr2str, x, digits, base)
}

formatMpfr <-
    function(x, digits = NULL, trim = FALSE, scientific = NA,
             base = 10, showNeg0 = TRUE,
	     big.mark = "", big.interval = 3L,
	     small.mark = "", small.interval = 5L,
             decimal.mark = ".",
             exponent.char = if(base <= 14) "e" else if(base <= 36) "E" else "|e",
	     zero.print = NULL, drop0trailing = FALSE, ...)
{
    ##	digits = NULL : use as many digits "as needed"
    ff <- .mpfr2str(x, digits, base=base)

    isNum <- ff$finite	## ff$finite == is.finite(x)
    i0 <- ff$is.0	## == mpfrIs0(x)
    ex <- ff$exp ## the *decimal* exp : one too large *unless* x == 0
    r  <- ff$str
    if(is.null(digits)) digits <- nchar(r)

    if(any(i0)) {
	## sign(x) == -1 "fails" for '-0'
	hasMinus <- substr(ff$str, 1L,1L) == "-"
	if(showNeg0) {
	} else if(any(iN0 <- hasMinus & i0)) {
	    ## get rid of "-" for "negative zero"
	    r[iN0] <- substring(r[iN0], 2)
	    hasMinus[iN0] <- FALSE
	}
	Ex <- ex
	Ex[!i0] <- ex[!i0] - 1L
    } else {
	Ex <- ex - 1L
	hasMinus <- sign(x) == -1
    }

    if(!all(isNum)) ## "@Inf@", "@NaN@", ...
	r[!isNum] <- gsub("@", '', r[!isNum], fixed=TRUE)

    ##' (maybe) add decimal point after position  k
    patch <- function(str, k)
	paste(substr   (str, 1L, k),
	      substring(str, k+1L), sep = decimal.mark)

    if(is.na(scientific))
	scientific <- as.numeric(getOption("scipen"))
    ## This very much depends on the desired format.
    ## if(scientific) --> all get a final "e<exp>"; otherwise, we
    ## adopt the following simple scheme :
    hasE <- { if(is.logical(scientific)) scientific else
	      isNum & (Ex < -4 + scientific | Ex > digits) }

    if(aE <- any(ii <- isNum & hasE)) {
        ii <- which(ii)
	i. <- 1L + hasMinus
	r[ii] <- patch(r[ii], i.[ii])
	if(drop0trailing)
	    ## drop 0's only *after* (and together with!) decimal mark:
	    r[ii] <- sub(paste0(decimal.mark, "0+$"), "", r[ii])
	r[ii] <- paste(r[ii], as.character(Ex[ii]), sep = exponent.char)
    }
    use.prettyN <- (base <= 14 && (!aE || exponent.char == "e"))
    if(non.sci <- !all(hasE)) { ## "non-scientific" i.e. without final  e<nn> :
	ii <- isNum & !hasE
	## iNeg <- ex <= 0 & ii ## i.e., ex	 in {0,-1,-2,-3}
	## iPos <- ex >  0 & ii ## i.e., ex	 in {1,2..., digits}
	iNeg <- Ex <  0	 & ii ## i.e., ex	 in {0,-1,-2,-3}
	iPos <- Ex >= 0	 & ii ## i.e., ex	 in {1,2..., digits}

	nZeros <- function(n) ## e.g.  nZeros(2:0) gives  c("00","0", "")
	    vapply(n, function(k) paste(rep.int("0", k), collapse = ""), "")
	if(any(eq <- (Ex == digits))) {
	    r[eq] <- paste0(r[eq], "0")
	    Ex[eq] <- Ex[eq] + 1L
	}
	if(any(iNeg)) { ## "0.00..." : be careful with minus sign
	    if(any(isMin <- hasMinus[iNeg])) {
		rr <- r[iNeg]
		rr[isMin] <- substring(rr[isMin], 2)
		r[iNeg] <- paste0(c("","-")[1+isMin], "0.",
				  nZeros(-ex[iNeg]), rr)
	    }
	    else {
		r[iNeg] <- paste0("0.", nZeros(-ex[iNeg]), r[iNeg])
	    }
	}
	if(any(iPos)) ## "xy.nnnn" :
	    r[iPos] <- patch(r[iPos], (hasMinus + Ex+1L)[iPos])
    }
    if(use.prettyN)
        r <- prettyNum(r, big.mark = big.mark, big.interval = big.interval,
                       small.mark = small.mark,
                       small.interval = small.interval,
                       decimal.mark = decimal.mark,
                       zero.print = zero.print, drop0trailing = drop0trailing,
                       preserve.width = if (trim) "individual" else "common")
    else {
	if(non.sci && drop0trailing)
	    ## drop 0's only *after* (and together with!) decimal mark:
	    r <- sub(paste0(decimal.mark, "0+$"), "", r)
	if(!missing(big.mark) || !missing(big.interval) || !missing(small.interval) ||
	    !missing(small.mark) || !missing(big.interval) || !missing(zero.print))
	    warning("with base >= 15 or 'exponent.char != \"e\", cannot use prettyNum()")
    }
    if(is.null(d <- dim(x))) r
    else array(r, dim=d, dimnames = dimnames(x))
}
setMethod("format", "mpfr", formatMpfr)

formatN.mpfr <- function(x, drop0trailing = TRUE, ...) {
    paste0(formatMpfr(x, drop0trailing=drop0trailing, ...),"_M")
}

setAs("mpfr", "character", function(from)
      format(from, digits=NULL, drop0trailing = TRUE))

setAs("character", "mpfr", function(from) mpfr(from))
