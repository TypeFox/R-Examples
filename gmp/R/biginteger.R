is.bigz <- function(x) is.raw(x) && inherits(x, "bigz")
is.bigq <- function(x) is.raw(x) && inherits(x, "bigq")

## Auxiliaries:
if(getRversion() < "2.15")
    paste0 <- function(...) paste(..., sep = '')

setGeneric("asNumeric", useAsDefault = function(x) {
    if(is.numeric(x)) x else if(is.atomic(x)) {
        storage.mode(x) <- "numeric"; x }
    else as(x, "numeric")
})

#----------------------------------------------------------
#
#  Author        : Immanuel Scholz (immanuel.scholz@gmx.de)
#		   Technische Universitaet Dresden
#
#  Brief         : Stub to call the dll functions
#
#  Licence       : GPL
#
#----------------------------------------------------------

"+.bigz" <- add.bigz <- function(e1, e2) .Call(biginteger_add, e1, e2)

"-.bigz" <- sub.bigz <- function(e1, e2=NULL)
{
    if(is.null(e2))
      .Call(biginteger_sub, 0, e1)
    else
      .Call(biginteger_sub, e1, e2)
}

"*.bigz" <- mul.bigz <- function(e1, e2) .Call(biginteger_mul, e1, e2)

## divq : integer division
"%/%.bigz" <- divq.bigz <- function(e1, e2) .Call(biginteger_divq, e1, e2)

## div : division of integers -> either rational or (mod) integer division
"/.bigz" <- div.bigz <- function(e1, e2) .Call(biginteger_div, e1, e2)

"%%.bigz" <- mod.bigz <- function(e1, e2) .Call(biginteger_mod,e1, e2)

"^.bigz" <- pow.bigz <- function(e1, e2,...) .Call(biginteger_pow,e1, e2)

inv.bigz <- function(a,b,...) .Call(biginteger_inv,a,b)

gcd <- function(a,b)
      UseMethod("gcd")
gcd.default <- function(a,b) as.integer(gcd.bigz(a,b))
gcd.bigz <- function(a,b) .Call(biginteger_gcd,a,b)

## just because lcm() is a trivial function in 'graphics' .. hmm
##lcm <- function(a,b)
##      UseMethod("lcm")

lcm.default <- function(a,b)
  as.integer(lcm.bigz(a,b))
lcm.bigz <- function(a,b) .Call(biginteger_lcm,a,b)

print.bigz <- function(x, quote = FALSE, initLine = is.null(modulus(x)), ...)
{
  if((n <- length(x)) > 0) {
    if(initLine) {
      cat("Big Integer ('bigz') ")
      kind <- if(isM <- !is.null(nr <- attr(x, "nrow")))
        sprintf("%d x %d matrix", nr, n/nr)
      else if(n > 1) sprintf("object of length %d", n) else ""
      cat(kind,":\n", sep="")
    }
    print(as.character(x), quote = quote, ...)
  }
  else
    cat("bigz(0)\n")
  invisible(x)
}

as.bigz <- function(a, mod = NA)
{
  if(is.null(mod)) mod <- NA
  if(inherits(a, "bigq"))
    as.bigz.bigq(a, mod)
  else
    .Call(biginteger_as, a, mod)
}
## the .as*() functions are exported for Rmpfr
.as.bigz <- function(a, mod = NA) {
  if(inherits(a, "bigq")) as.bigz.bigq(a, mod) else .Call(biginteger_as, a, mod)
}
..as.bigz <- function(a, mod = NA) .Call(biginteger_as, a, mod)

.as.char.bigz <-
as.character.bigz <-
    function(x, b = 10L, ...) .Call(biginteger_as_character, x, b)


##' format() Numbers such as to distinguish  bigz, integer, double, mpfr, etc
formatN <- function(x, ...) UseMethod("formatN")
formatN.integer <- function(x, ...) paste0(as.character(x, ...), "L")
formatN.bigz    <- function(x, ...) {
    r <- as.character(x, ...)
    if(any(noMod <- is.null(modulus(x))))
	r[noMod] <- paste0(r[noMod],"_Z")
    r
}
formatN.double	<- function(x, ...) {
    r <- vapply(x, format, "", ...)
    if(any(intLike <- !grepl("[^-0-9]",r)))
	r[intLike] <- paste0(r[intLike],".")
    r
}
##' Default Method: Use the standard format() --- e.g. for complex
formatN.default <- function(x, ...) format(x, ...)


as.double.bigz  <- function(x,...) .Call(biginteger_as_numeric, x)
as.integer.bigz <- function(x,...) .Call(biginteger_as_integer, x)

.bigz2num <- function(x) {
    r <- .Call(biginteger_as_numeric, x)
    if(!is.null(d <- dim(x))) dim(r) <- d
    r
}
setMethod("asNumeric", "bigz", .bigz2num)


length.bigz <- function(x) .Call(biginteger_length, x)

"length<-.bigz"<- function(x, value) .Call(biginteger_setlength, x, value)

modulus      <- function(a) UseMethod("modulus")
modulus.bigz <- function(a) attr(a, "mod")

`modulus<-`      <- function(a, value) UseMethod("modulus<-")
`modulus<-.bigz` <- function(a, value) as.bigz(a, value)


## inv <- function(a,...) UseMethod("inv")

## pow <- function(a,...) UseMethod("pow")

powm <- function(x,y, n) .Call(biginteger_powm, x,y,n)

"<.bigz"  <- function(e1, e2) .Call(biginteger_lt, e1, e2)
">.bigz"  <- function(e1, e2) .Call(biginteger_gt, e1, e2)
"<=.bigz" <- function(e1, e2) .Call(biginteger_lte, e1, e2)
">=.bigz" <- function(e1, e2) .Call(biginteger_gte, e1, e2)
"==.bigz" <- function(e1, e2) .Call(biginteger_eq, e1, e2)
"!=.bigz" <- function(e1, e2) .Call(biginteger_neq, e1, e2)

is.whole <- function(x) UseMethod("is.whole")
is.whole.default <- function(x) {
    n <- length(x)
    if(is.atomic(x)) {
        if(is.integer(x) || is.logical(x))
            return(rep.int(TRUE, n))
        if(is.numeric(x))
            return(x == floor(x))
        if(is.complex(x))
            return(x == round(x))
    }
    ## else:
    logical(n) ## == rep.int(FALSE, length(x))
}


is.na.bigz <- function(x) .Call(biginteger_is_na, x)
is.finite.bigz <- function(x) !is.na.bigz(x) # otherwise all are finite
is.whole.bigz <- function(x) rep.int(TRUE, length(x))
is.infinite.bigz <- function(x) rep.int(FALSE, length(x))

frexpZ <- function(x) .Call(bigI_frexp, x)

##' @title log2(Inverse of frexpZ(a))
##' @param L list(d = ., exp = .)
##' @return numeric vector
##' @author Martin Maechler
lg2.invFrexp <- function(L) {
    stopifnot(is.list(L), is.numeric(d <- L$d), is.numeric(ex <- L$exp),
	      (n <- length(d)) == length(ex))
    ex + log2(d)
}

###------------------------- 'Math' S3 group ------------------------------

## Most 'Math' group would be hard to implement --- [TODO via Rmpfr -- or stop("...via Rmpfr")?
## Fall-back: *not* implemented  {or use as.double() ??}
Math.bigz <- function(x, ...) { .NotYetImplemented() }
            ## • ‘abs’, ‘sign’, ‘sqrt’,
            ##   ‘floor’, ‘ceiling’, ‘trunc’,
            ##   ‘round’, ‘signif’
            ## • ‘exp’, ‘log’, ‘expm1’, ‘log1p’,
            ##   ‘cos’, ‘sin’, ‘tan’,
            ##   ‘acos’, ‘asin’, ‘atan’
            ##   ‘cosh’, ‘sinh’, ‘tanh’,
            ##   ‘acosh’, ‘asinh’, ‘atanh’
            ## • ‘lgamma’, ‘gamma’, ‘digamma’, ‘trigamma’
            ## • ‘cumsum’, ‘cumprod’, ‘cummax’, ‘cummin’

abs.bigz <- function(x) .Call(biginteger_abs,x)
sign.bigz <- function(x) .Call(biginteger_sgn,x)

floor.bigz <- ceiling.bigz <- function(x) x
trunc.bigz <- function(x, ...) x
round.bigz <- function(x, digits=0) {
    if(digits == 0) x
    else stop("digits != 0  is not yet implemented")
}

gamma.bigz <- function(x) factorialZ(x-1)

cumsum.bigz <- function(x) .Call(biginteger_cumsum, x)
## TODO: add cummax(), cummin(), cumprod()


log2.bigz <- function(x) .Call(biginteger_log2, x)
## not exported:
ln.bigz	 <- function(x) .Call(biginteger_log, x)
log.bigz <- function(x, base=exp(1))
{
    if(missing(base))
	ln.bigz(x)
    else
	ln.bigz(x)/log(base)
}

log10.bigz <- function(x) ln.bigz(x) / log(10)

##------------end{'Math'} group -------------------------------------


###------------------------- 'Summary' S3 group ------------------------------
##---- "max"   "min"   "range" "prod"  "sum"   "any"   "all"  -----

max.bigz <- function(...,na.rm=FALSE)
{
 .Call(biginteger_max, c.bigz(...), na.rm)
}

min.bigz <- function(...,na.rm=FALSE)
{
 .Call(biginteger_min, c.bigz(...), na.rm)
}

## range(): works automatically via  range.default() and the above min(), max()

prod.bigz <- function(..., na.rm = FALSE)
{
    X <- c.bigz(...)
   .Call(biginteger_prod, if(na.rm) X[!is.na(X)] else X)
}

sum.bigz <- function(..., na.rm = FALSE)
{
    X <- c.bigz(...)
    .Call(biginteger_sum, if(na.rm) X[!is.na(X)] else X)
}
##------------end{Summary group}------------------------------------

## FIXME: implement faster in C
setMethod("which.max", "bigz", function(x) which.max(x == max(x)))
setMethod("which.min", "bigz", function(x) which.max(x == min(x)))

c.bigz <- function(..., recursive = FALSE)
{
    .Call(biginteger_c, list(...))
}

## This is practically identical to  grid :: rep.unit :
rep.bigz <- function(x, times=1, length.out=NA, each=1, ...) {
    ## if (length(x) == 0)
    ##	   stop("invalid 'unit' object")
    if(!missing(times) && missing(length.out) && missing(each))
        .Call(biginteger_rep, x, times)
    else {
	## Determine an appropriate index, then call subsetting code
	x[ rep(seq_along(x), times=times, length.out=length.out, each=each) ]
    }
}

duplicated.bigz <- function(x, incomparables = FALSE, ...) {
    x <- as.character(x) # lazy and inefficient --> TODO in C++
    NextMethod("duplicated", x)
}

unique.bigz <- function(x, incomparables = FALSE, ...)
    x[!duplicated(x, incomparables, ...)]

all.equal.bigz <- function(target, current, ...)
    if(all(target == current)) TRUE else "'target'(bigz) and 'current' differ"

# Isprime, return:
#   0 if not prime
#   1 if probably prime
#   2 if prime
isprime <- function(n,reps=40)
  {
    .Call(biginteger_is_prime, n, as.integer(reps))
  }

nextprime <- function(n) .Call(biginteger_nextprime, n)

gcdex <- function(a, b) .Call(biginteger_gcdex, a, b)

urand.bigz <- function(nb=1, size=200, seed=0)
  {
    ok <- (seed != 0)
    .Call(biginteger_rand_u,
          as.integer(nb),
          as.integer(size),
          seed,
          as.integer(ok))
  }

sizeinbase <- function(a, b=10)
{
    if(as.integer(b) < 2)
        stop("base must be >= 2")
    .Call(biginteger_sizeinbase, a, as.integer(b))
}

factorialZ  <- function(n) .Call(bigI_factorial, as.integer(n))
chooseZ  <- function(n, k) .Call(bigI_choose,  n, as.integer(k))

fibnum  <- function(n) .Call(bigI_fibnum,  as.integer(n))
fibnum2 <- function(n) .Call(bigI_fibnum2, as.integer(n))

lucnum  <- function(n) .Call(bigI_lucnum,  as.integer(n))
lucnum2 <- function(n) .Call(bigI_lucnum2, as.integer(n))

factorize <- function(n) .Call(factorR, as.bigz(n))

solve.bigz <- function(a, b,...)
  {
    if(missing(b))
      .Call(inverse_z, a)
    else
      .Call(solve_z, a, b)
  }


`[[.bigz` <- function(x, i=NA) .Call(biginteger_get_at, x, i)

`[[<-.bigz` <- function(x, i=NA, value)
    .Call(biginteger_set_at, x, i, value)

`[.bigz` <- function(x, i=NULL, j=NULL, drop=TRUE)
{
  mdrop <- missing(drop)
  Narg <- nargs() - (!mdrop)
  has.j <- !missing(j)
  if(isM <- !is.null(nr <- attr(x, "nrow"))) { ## matrix
    ## FIXME  x[i,] vs. x[,j] vs. x[i]
    .Call(matrix_get_at_z, x, i,j)
  } else { ## non-matrix
    if(has.j) stop("invalid vector subsetting")
    ## ugly "workaround"
    r <- .Call(matrix_get_at_z, x, i, NULL)
    attr(r,"nrow") <- NULL
    r
  }

}

`[<-.bigz` <- function(x, i=NULL, j=NULL, value)
{
  .Call(matrix_set_at_z, x, value, i,j)
}


