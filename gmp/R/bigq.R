#----------------------------------------------------------
#
#  Author        : Antoine Lucas (adapted from biginteger class made by
#                                 Immanuel Scholz)
#
#  Brief         : Stub to call the dll functions
#
#  Licence       : GPL
#
#----------------------------------------------------------


add.bigq <- `+.bigq` <- function(e1, e2) .Call(bigrational_add, e1, e2)

sub.bigq<- `-.bigq` <- function(e1, e2=NULL) {
  if(is.null(e2))
    .Call(bigrational_sub,0,e1)
  else
    .Call(bigrational_sub,e1,e2)
}

mul.bigq <- `*.bigq` <- function(e1, e2) .Call(bigrational_mul, e1, e2)

"/.bigq" <- div.bigq <- function(e1, e2) .Call(bigrational_div, e1, e2)

## note: Here, the 2nd argument must be big integer, not rational:
"^.bigq" <- pow.bigq <- function(e1, e2) .Call(bigrational_pow, e1, as.bigz(e2))

print.bigq <- function(x, quote = FALSE, initLine = TRUE, ...)
{
  if((n <- length(x)) > 0) {
    if(initLine) {
      cat("Big Rational ('bigq') ")
      kind <- if(isM <- !is.null(nr <- attr(x, "nrow")))
        sprintf("%d x %d matrix", nr, n/nr)
      else if(n > 1) sprintf("object of length %d", n) else ""
      cat(kind,":\n", sep="")
    }
    print(as.character(x), quote = quote, ...)
  }
  else
    cat("bigq(0)\n")
  invisible(x)
}

as.bigq <- function(n, d=1)
{
    .Call(bigrational_as, n, d)
}

as.character.bigq <- function(x, b = 10L, ...)
{
    .Call(bigrational_as_character, x, b)
}

formatN.bigq	<- function(x, ...) {
    r <- as.character(x, ...)
    if(any(iI <- is.whole.bigq(x)))
	r[iI] <- paste0(r[iI],"/1")
    r
}

as.double.bigq<- function(x,...) .Call(bigrational_as_numeric, x)
## maybe sub-optimal, but at least "R-consistent" in warnings/errors...:
as.integer.bigq<- function(x,...) as.integer(.Call(bigrational_as_numeric, x))

.bigq2num <- function(x) {
    ## cat(".bigq2num():\n")
    r <- .Call(bigrational_as_numeric, x)
    if(!is.null(d <- dim(x))) dim(r) <- d
    r
}
setMethod("asNumeric", "bigq", .bigq2num)


denominator <- function(x) {
  r <- .Call(bigrational_den,x)
  if(!is.null(d <- dim(x))) dim(r) <- d
  r
}

"denominator<-" <- function(x,value)
  as.bigq(numerator(x),value)

numerator <- function(x) {
    r <- .Call(bigrational_num,x)
    if(!is.null(d <- dim(x))) dim(r) <- d
    r
}

"numerator<-" <- function(x,value)
  as.bigq(value,denominator(x))


as.bigz.bigq <- function(a, mod = NA)
{
  ## "FIXME":  considerably faster in C / C++
  if(any(ina <- is.na.bigq(a))) {
    r <- as.bigz(rep.int(NA, length(a)))
    if(any(ii <- !ina)) {
	a <- a[ii]
	r[ii] <- as.bigz(numerator(a) %/% denominator(a), mod[ii])
    }
    attr(r,"nrow") <- attr(a, "nrow")
    r
  }
  else # no NA's
    as.bigz(numerator(a) %/% denominator(a), mod)
}

length.bigq<- function(x) .Call(bigrational_length, x)
`length<-.bigq` <- function(x, value) .Call(bigrational_setlength, x, value)


"<.bigq"  <- function(e1,e2) .Call(bigrational_lt,  e1, e2)
">.bigq"  <- function(e1,e2) .Call(bigrational_gt,  e1, e2)
"<=.bigq" <- function(e1,e2) .Call(bigrational_lte, e1, e2)
">=.bigq" <- function(e1,e2) .Call(bigrational_gte, e1, e2)
"==.bigq" <- function(e1,e2) .Call(bigrational_eq,  e1, e2)
"!=.bigq" <- function(e1,e2) .Call(bigrational_neq, e1, e2)

is.na.bigq <- function(x) .Call(bigrational_is_na, x)
is.whole.bigq <- function(x) .Call(bigrational_is_int, x)
is.finite.bigq <- function(x) !is.na.bigq(x) # otherwise all are finite
is.infinite.bigq <- function(x) rep.int(FALSE, length(x))

if(FALSE) ## This does not work: is.atomic is primitive *NON*-generic:
is.atomic.bigq <- function(x) FALSE # otherwise does return TRUE !

###  <bigz> o <bigq>  --- really dispatch on two arguments --> use S4
if(FALSE) { ## not working really --- see also ./matrix-prods.R

setMethod("Ops", signature(e1 = "bigq", e2 = "bigz"),
	  function(e1, e2) callGeneric(e1, as.bigq(e2)))
setMethod("Ops", signature(e1 = "bigz", e2 = "bigq"),
	  function(e1, e2) callGeneric(as.bigq(e1), e2))
}

###------------------------- 'Math' S3 group ------------------------------

## Most 'Math' group would be hard to implement --- [TODO via Rmpfr -- or stop("...via Rmpfr")?
## Fall-back: *not* implemented  {or use as.double() ??}
Math.bigq <- function(x, ...) { .NotYetImplemented() }


abs.bigq <- function(x) {
    numerator(x) <- abs(numerator(x))
    x
}

sign.bigq <- function(x) sign(numerator(x))

## TODO
trunc.bigq <- function(x, ...) .NotYetImplemented()
floor.bigq   <- function(x) .NotYetImplemented()
ceiling.bigq <- function(x) .NotYetImplemented()
round.bigq <- function(x, digits = 0) {
    .NotYetImplemented()
}

cumsum.bigq <- function(x) .Call(bigrational_cumsum, x)
## TODO: add cummax(), cummin(), cumprod()

## FIXME: implement  log() etc --- see ./biginteger.R

##------------end{'Math'} group -------------------------------------




c.bigq <- function(..., recursive = FALSE) {
    .Call(bigrational_c, list(...))
}

## This is practically identical to  grid :: rep.unit :
rep.bigq <- function(x, times=1, length.out=NA, each=1, ...) {
    ## if (length(x) == 0)
    ##	   stop("invalid 'unit' object")
    if(!missing(times) && missing(length.out) && missing(each))
	.Call(bigrational_rep,x,times)
    else {
	## Determine an appropriate index, then call subsetting code
	x[ rep(seq_along(x), times=times, length.out=length.out, each=each) ]
    }
}

duplicated.bigq <- duplicated.bigz ## cheap (for now)

## unique.bigq <- function(x, incomparables = FALSE, ...)
##     x[!duplicated(x, incomparables=incomparables, ...)]
unique.bigq <- unique.bigz

##' mean() method needed for all.equal.bigq() below:
mean.bigq <- function(x, trim = 0, na.rm = FALSE, ...) {
    if(trim != 0) stop("'trim > 0' is not yet implemented for \"bigq\"")
    if (na.rm) x <- x[!is.na(x)]
    sum(x) / length(x)
}

## Almost:
## all.equal.bigq <- all.equal.numeric
## environment(all.equal.bigq) <- environment()# i.e. of 'gmp'  name space
## but we copy-paste all.equal.numeric  {and slightly modify}:
all.equal.bigq <-
    function(target, current, tolerance = .Machine$double.eps ^ .5,
             scale = NULL, check.attributes = FALSE, check.class=FALSE, ...)
{
    msg <- if(check.attributes)
	attr.all.equal(target, current, tolerance=tolerance, scale=scale, ...)
    if(check.class && data.class(target) != data.class(current)) {
	msg <- c(msg, paste0("target is ", data.class(target), ", current is ",
                             data.class(current)))
	return(msg)
    }

    lt <- length(target)
    lc <- length(current)
    cplx <- FALSE
    if(lt != lc) {
	## *replace* the 'Lengths' msg[] from attr.all.equal():
	if(!is.null(msg)) msg <- msg[- grep("\\bLengths\\b", msg)]
	msg <- c(msg, paste0("bigq",
                             ": lengths (", lt, ", ", lc, ") differ"))
	return(msg)
    }
    ## remove atttributes (remember these are both numeric or complex vectors)
    ## one place this is needed is to unclass Surv objects in the rpart test suite.
    target <- as.vector(target)
    current <- as.vector(current)
    out <- is.na(target)
    if(any(out != is.na(current))) {
	msg <- c(msg, paste("'is.NA' value mismatch:", sum(is.na(current)),
			    "in current", sum(out), "in target"))
	return(msg)
    }
    out <- out | target == current
    if(all(out)) { if (is.null(msg)) return(TRUE) else return(msg) }

    target <- target[!out]
    current <- current[!out]
    if(is.integer(target) && is.integer(current)) target <- as.double(target)
    xy <- mean((if(cplx) Mod else abs)(target - current))
    what <-
	if(is.null(scale)) {
	    xn <- mean(abs(target))
	    if(is.finite(xn) && xn > tolerance) {
		xy <- xy/xn
		"relative"
	    } else "absolute"
	} else {
	    xy <- xy/scale
	    "scaled"
	}

    if (cplx) what <- paste(what, "Mod") # PR#10575
    if(is.na(xy) || xy > tolerance)
        msg <- c(msg, paste("Mean", what, "difference:", format(xy)))

    if(is.null(msg)) TRUE else msg
}



solve.bigq <- function(a,b,...)
  {
    if(missing(b))
      .Call(inverse_q,a)
    else
      .Call(solve_q,a,b)
  }


`[[.bigq`<- function(x, i=NA)
{
    .Call(bigrational_get_at, x, i)
}

`[[<-.bigq` <- function(x, i=NA, value)
{
    .Call(bigrational_set_at, x, i, value)
}


`[.bigq` <- function(x,i=NULL,j=NULL, drop=TRUE)
{
  .Call(matrix_get_at_q, x, i,j)
}


`[<-.bigq` <- function(x,i=NULL,j=NULL,value)
{
  .Call(matrix_set_at_q, x, value,i,j )
}


max.bigq <- function(...,na.rm=FALSE)
{
 .Call(bigrational_max, c.bigq(...), na.rm)
}

min.bigq <- function(...,na.rm=FALSE)
{
 .Call(bigrational_min, c.bigq(...), na.rm)
}

## FIXME: implement faster in C
setMethod("which.max", "bigq", function(x) which.max(x == max(x)))
setMethod("which.min", "bigq", function(x) which.max(x == min(x)))

sum.bigq <- function(..., na.rm = FALSE)
{
    X <- c.bigq(...)
   .Call(bigrational_sum, if(na.rm) X[!is.na(X)] else X)
}

prod.bigq <- function(..., na.rm = FALSE)
{
    X <- c.bigq(...)
   .Call(bigrational_prod, if(na.rm) X[!is.na(X)] else X)
}



