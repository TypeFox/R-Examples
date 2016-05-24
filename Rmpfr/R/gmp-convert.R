#### Conversions   bigz  <-> mpfr   // also bigq <--> mpfr

if(packageVersion("gmp") < "0.5.8")## <-> ../NAMESPACE
    is.matrixZQ <- function(x) !is.null(attr(x, "nrow"))

## The following code is experimental, hence the "." :

### FIXME: we go via character.. which is not really efficient.
## ------  rather "should" use  MPFR Functions
## int mpfr_set_z (mpfr_t ROP, mpz_t OP, mpfr_rnd_t RND)
## int mpfr_set_q (mpfr_t ROP, mpq_t OP, mpfr_rnd_t RND)
##
## Set the value of ROP from OP, rounded toward the given direction RND.
##
## Directly in C, we'd need both Rmpfr and gmp's  C code (!)
## TODO(?:  gmp should "export" its C++ API ( -> inst/include/*.hh )
## and we should add  'LinkingTo: gmp' to DESCRIPTION and
##  then use C++ with "C" { ...} for those parts
.bigz2mpfr <- function(x, precB = NULL, rnd.mode = c('N','D','U','Z','A')) {
    stopifnot(inherits(x, "bigz"))
    ..bigz2mpfr(x, precB, rnd.mode)
}
## Fast, no-checking (and not exported) version:
..bigz2mpfr <- function(x, precB = NULL, rnd.mode = c('N','D','U','Z','A'))
    ## precB: 4 == log2(16) = log(base)
{
    stopifnot(is.character(rnd.mode <- toupper(rnd.mode)))
    rnd.mode <- match.arg(rnd.mode)
    b <- 16L
    cx <- .as.char.bigz(x, b)
    if(is.null(precB)) precB <- 4L*nchar(cx)
    if(is.matrixZQ(x))
	new("mpfrMatrix", .Call(str2mpfr1_list, cx, precB, b, rnd.mode),
	    Dim = dim(x))# "bigz" has no dimnames
    else
	new("mpfr", .Call(str2mpfr1_list, cx, precB, b, rnd.mode))
}
setAs("bigz", "mpfr", function(from) ..bigz2mpfr(from))


## FIXME: rather should use MPFR -- Function :
## ----   int mpfr_get_z (mpz_t ROP, mpfr_t OP, mpfr_rnd_t RND)
## Convert OP to a `mpz_t', after rounding it with respect to RND.  ....
## FIXME(2): should 'gmp' change as.bigz into an S3 generic, so this becomes S3 method?
as.bigz.mpfr <-
.mpfr2bigz <- function(x, mod=NA) {
    if(is.null(mod)) mod <- NA_integer_
    stopifnot(is(x, "mpfr"),
	      is.na(mod) || (length(mod) == 1L && is.numeric(mod)))
    dx <- dim(x)
### FIXME or rather  roundMpfr()  [or even round "RND" as in mpfr_get_z() above] ??
    cx <- format(trunc(x), drop0trailing=TRUE)
    dim(cx) <- dx ## needed?? {should *not* be, as in base R!}
    ..as.bigz(cx, mod)
}
setAs("mpfr", "bigz", function(from) .mpfr2bigz(from))


## Fast, no-checking (and not exported) version:
..bigq2mpfr <- function(x, precB = NULL, rnd.mode = c('N','D','U','Z','A')) {
    stopifnot(is.character(rnd.mode <- toupper(rnd.mode)))
    rnd.mode <- match.arg(rnd.mode)
    N <- numerator(x)
    D <- denominator(x)
    if(is.null(precB)) {
        eN <- frexpZ(N)$exp
        eD <- frexpZ(D)$exp
        precB <- pmax(128L, eN + eD + 1L) # precision of result
    }
    ..bigz2mpfr(N, precB, rnd.mode) / ..bigz2mpfr(D, precB, rnd.mode)
}

.bigq2mpfr <- function(x, precB = NULL, rnd.mode = c('N','D','U','Z','A')) {
    stopifnot(inherits(x, "bigq"))
    ..bigq2mpfr(x, precB, rnd.mode)
}
setAs("bigq", "mpfr", function(from) ..bigq2mpfr(from))

## TODO(?)  "mpfr" ->  "bigq"
## a) in the spirit  MASS::fractions()    or
## b) "native" MPFR "support" -- not yet available: has mpfr_get_z() but not get_q()
