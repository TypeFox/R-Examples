#### Define mpfr methods for Math  and Math2  group functions
####                        ======     =====

### "Arith", "Compare",..., are in ./Arith.R
###  ----                            ~~~~~~~

## [1] "abs"    "sign"    "sqrt"    "ceiling" "floor" "trunc" "cummax"
## [8] "cummin" "cumprod" "cumsum"  "exp"     "expm1" "log"   "log10"
##[15] "log2"   "log1p"   "cos"     "cosh"    "sin"   "sinh"  "tan"
##[22] "tanh"   "acos"    "acosh"   "asin"    "asinh" "atan"  "atanh"
##[29] "gamma"  "lgamma"  "digamma" "trigamma"

if(FALSE) ## here are the individual function
    dput(getGroupMembers("Math"))

## Uniform interface to C:
##
## Pass integer code to call and do the rest in C
## Codes from ~/R/D/r-devel/R/src/main/names.c :
.Math.codes <-
    c(
      "floor" =     1,
      "ceiling" =   2,
      "sqrt" =      3,
      "sign" =      4,
      "exp" =      10,
      "expm1" =    11,
      "log1p" =    12,
      "cos" =      20,
      "sin" =      21,
      "tan" =      22,
      "acos" =     23,
      "asin" =     24,
      "cosh" =     30,
      "sinh" =     31,
      "tanh" =     32,
      "acosh" =    33,
      "asinh" =    34,
      "atanh" =    35,
      "lgamma" =   40,
      "gamma" =    41,
      "digamma" =  42,
      "trigamma" = 43,
	## R >= 3.1.0 :
      "cospi" =    47,
      "sinpi" =    48,
      "tanpi" =    49
        )

.Math.gen <- getGroupMembers("Math")

## Those "Math" group generics that are not in the do_math1 table above

.Math.codes <-
    c(.Math.codes,
      "trunc" = 0, "atan" = 25, # "abs" has own method!
      "log" = 13, "log2" = 14, "log10" = 15,
      "cummax" = 71, "cummin" = 72, "cumprod" = 73, "cumsum" = 74,
      ## These are *NOT* in R's  Math group, but 1-argument math functions
      ## available in the mpfr - library:
      "erf" = 101, "erfc" = 102, "zeta" = 104, "Eint" = 106, "Li2" = 107,
      "j0" = 111, "j1" = 112, "y0" = 113, "y1" = 114,
      "Ai" = 120) # Airy function (new in mpfr 3.0.0)
storage.mode(.Math.codes) <- "integer"

if(FALSE)
.Math.gen[!(.Math.gen %in% names(.Math.codes))]
## "abs" -- only one left

## A few ones have a very simple method:
## Note that the 'sign' slot is from the C-internal struct
## and is always +/- 1 , but R's sign(0) |--> 0
.getSign <- function(x) vapply(getD(x), slot, 1L, "sign")
.mpfr.sign <- function(x) {
    r <- numeric(length(x))# all 0
    not0 <- !mpfrIs0(x)
    r[not0] <- .getSign(x[not0])
    r
}
setMethod("sign", "mpfr", .mpfr.sign)

## R version, no longer used:
.abs.mpfr <- function(x) {
    ## FIXME: faster if this happened in a .Call
    xD <- getDataPart(x)   # << currently [2011] *faster* than  x@Data
    for(i in seq_along(x))
        slot(xD[[i]], "sign", check=FALSE) <- 1L
    setDataPart(x, xD, check=FALSE) ## faster than  x@.Data <- xD
}
setMethod("abs", "mpfr",
	  function(x) .Call(Rmpfr_abs, x))

## Simple methods for "complex" numbers, just so "they work"
setMethod("Re",  "mpfr", function(z) z)
setMethod("Im",  "mpfr", function(z) 0*z)
setMethod("Conj","mpfr", function(z) z)
setMethod("Mod", "mpfr", function(z) abs(z))
setMethod("Arg", "mpfr", function(z) {
    prec <- .getPrec(z)
    r <- mpfr(0, prec)
    neg <- !mpfrIs0(z) & .getSign(z) == -1
    r[neg] <- Const("pi", prec = prec[neg])
    r
})



## Note that  factorial() and lfactorial() automagically work through  [l]gamma()
## but for the sake of "exact for integer"
setMethod("factorial", "mpfr",
	  function(x) {
	      r <- gamma(x + 1)
	      isi <- .mpfr.is.whole(x)
	      r[isi] <- round(r[isi])
	      r
	  })
## The "real" thing is to use  the MPFR-internal function:
factorialMpfr <- function(n, precBits = max(2, ceiling(lgamma(n+1)/log(2))),
                          rnd.mode = c('N','D','U','Z','A'))
{
    stopifnot(n >= 0)
    new("mpfr", .Call(R_mpfr_fac, n, precBits, match.arg(rnd.mode)))
}

##' Pochhammer rising factorial = Pochhammer(a,n) {1 of 2 definitions!}
##' we use the *rising* factorial for Pochhamer(a,n), i.e.,
##' the definition that the GSL and Mathematica use as well.
##' We want to do this well for *integer* n, only the general case is using
##' P(a,x) := Gamma(a+x)/Gamma(x)
pochMpfr <- function(a, n, rnd.mode = c('N','D','U','Z','A')) {
    stopifnot(n >= 0)
    if(!is(a, "mpfr")) ## use a high enough default precision (and recycle ..)
        a <- mpfr(a, precBits = pmax(1,n)*getPrec(a))
    else if((ln <- length(n)) != 1 && ln != length(a))
	a <- a + 0*n
    ## a@.Data[] <- .Call(R_mpfr_poch, a, n)
    ## a
    setDataPart(a, .Call(R_mpfr_poch, a, n, match.arg(rnd.mode)))
}

##' Binomial Coefficient choose(a,n)
##' We want to do this well for *integer* n
chooseMpfr <- function(a, n, rnd.mode = c('N','D','U','Z','A')) {
    stopifnot(n >= 0)
    if(!is(a, "mpfr")) { ## use high enough default precision
        lc <- lchoose(a,n)
        precB <- if(any(iF <- is.finite(lc))) ceiling(max(lc[iF])/log(2)) else 0
        ## add n bits for the n multiplications (and recycle {a,n} to same length)
        a <- mpfr(a, precBits = n + max(2, precB))
    } else if((ln <- length(n)) != 1 && ln != length(a))
	a <- a + 0*n
    ## a@.Data[] <- .Call(R_mpfr_choose, a, n)
    ## a
    setDataPart(a, .Call(R_mpfr_choose, a, n, match.arg(rnd.mode)))
}

chooseMpfr.all <- function(n, precBits=NULL, k0=1, alternating=FALSE) {
    ## return   chooseMpfr(n, k0:n) or (-1)^k * choose...  "but smartly"
    if(!is.numeric(n) || (n <- as.integer(n)) < 1)
	stop("n must be integer >= 1")
    stopifnot(is.numeric(n. <- k0), n. == (k0 <- as.integer(k0)),
              k0 <= n)
    sig <- if(alternating) (-1)^(k0:n) else rep.int(1, (n-k0+1))
    if(n == 1) return(mpfr(sig, 32))
    ## else : n >= 2
    n2 <- n %/% 2 # >= 1
    prec <- ceiling(lchoose(n,n2)/log(2)) # number of bits needed in result
    precBxtr <- max(2, n2 + prec) # need more for cumprod(), and division
    n2. <- mpfr(n2, precBxtr)
    r <- cumprod(seqMpfr(mpfr(n, precBxtr), n+1-n2., length.out=n2)) /
        cumprod(seqMpfr(1, n2., length.out=n2))
    prec <- max(2,prec)
    if(is.numeric(precBits) && (pB <- as.integer(round(precBits))) > prec)
	prec <- pB
    r <- roundMpfr(r, precBits = prec)
    ##
    ii <- c(seq_len(n2-1+(n%%2)), n2:1)
    if(k0 >= 2) ii <- ii[-seq_len(k0 - 1)]
    one <- .d2mpfr1(1, precBits=prec)
    r <- c(if(k0 == 0) one, getD(r)[ii], one)
    if(alternating) {
	for(i in seq_along(r)) if(sig[i] == -1)
	    slot(r[[i]], "sign", check=FALSE) <- - 1L
    }
    new("mpfr", r)
}## {chooseMpfr.all}


## http://en.wikipedia.org/wiki/N%C3%B6rlund%E2%80%93Rice_integral
## also deals with these alternating binomial sums
##'
##' version 1: already using the 'alternating' arg in chooseMpfr.all()
sumBinomMpfr.v1 <- function(n, f, n0=0, alternating=TRUE, precBits = 256)
{
    ## Note: n0 = 0, or 1 is typical, and hence chooseMpfr.all() makes sense
    stopifnot(0 <= n0, n0 <= n, is.function(f))
    sum(chooseMpfr.all(n, k0=n0, alternating=alternating) *
        f(mpfr(n0:n, precBits=precBits)))
}
##' version 2: chooseZ()*(-1)^(.) is considerably faster than chooseMpfr.all()
sumBinomMpfr.v2 <- function(n, f, n0=0, alternating=TRUE, precBits = 256,
			 f.k = f(mpfr(k, precBits=precBits)))
{
    ## Note: n0 = 0, or 1 is typical..
    stopifnot(0 <= n0, n0 <= n,
	      is.function(f) || (is(f.k, "mpfr") && length(f.k) == n-n0+1))
    k <- n0:n
    sum(if(alternating) chooseZ(n, k) * (-1)^(n-k) * f.k
        else chooseZ(n, k) * f.k)
}
## NB:  pbetaI() in  ./special-fun.R  uses a special version..
## --- if we do this *fast* in C --  do  pbetaI()  as well.
sumBinomMpfr <- sumBinomMpfr.v2


##' Rounding to binary bits, not decimal digits. Closer to the number
##' representation, this also allows to increase or decrease a number's precBits
##' @title Rounding to binary bits, "mpfr-internally"
##' @param x an mpfr number (vector)
##' @param precBits integer specifying the desired precision in bits.
##' @return an mpfr number as \code{x} but with the new 'precBits' precision
##' @author Martin Maechler
roundMpfr <- function(x, precBits, rnd.mode = c('N','D','U','Z','A')) {
    stopifnot(is(x, "mpfr"))
    setDataPart(x, .Call(R_mpfr_round, x, precBits, match.arg(rnd.mode)))
}

## "log" is still special with its 'base' :
setMethod("log", signature(x = "mpfr"),
	  function(x, base) {
	      if(!missing(base) && base != exp(1))
		  stop("base != exp(1) is not yet implemented")
	      setDataPart(x, .Call(Math_mpfr, x, .Math.codes[["log"]]))
	  })

setMethod("Math", signature(x = "mpfr"), function(x)
	  setDataPart(x, .Call(Math_mpfr, x, .Math.codes[[.Generic]])))

setMethod("Math2", signature(x = "mpfr"),
	  function(x, digits) {
	      ## NOTA BENE: vectorized in  'x'
	      if(any(ret.x <- !is.finite(x) | mpfrIs0(x))) {
		  if(any(ok <- !ret.x))
		      x[ok] <- callGeneric(x[ok], digits=digits)
		  return(x)
	      }
              if(!missing(digits)) {
                  digits <- as.integer(round(digits))
                  if(is.na(digits)) return(x + digits)
              } ## else: default *depends* on the generic

	      ## now: both x and digits are finite
	      pow10 <- function(d) mpfr(rep.int(10., length(d)),
					precBits = ceiling(log2(10)*as.numeric(d)))^ d
	      rint <- function(x) { ## have x >= 0 here
		  sml.x <- (x < .Machine$integer.max)
		  r <- x
		  if(any(sml.x)) {
		      x.5 <- x[sml.x] + 0.5
		      ix <- as.integer(x.5)
		      ## implement "round to even" :
		      if(any(doDec <- (abs(x.5 - ix) < 10*.Machine$double.eps & (ix %% 2))))
			  ix[doDec] <- ix[doDec] - 1L
		      r[sml.x] <- ix
		  }
		  if(!all(sml.x)) { ## large x - no longer care for round to even
		      r[!sml.x] <- floor(x[!sml.x] + 0.5)
		  }
		  r
	      }
	      neg.x <- x < 0
	      x[neg.x] <- - x[neg.x]
	      sgn <- ifelse(neg.x, -1, +1)
	      switch(.Generic,
		     "round" = { ## following ~/R/D/r-devel/R/src/nmath/fround.c :
			 if(missing(digits) || digits == 0)
			     sgn * rint(x)
			 else if(digits > 0) {
			     p10 <- pow10(digits)
			     intx <- floor(x)
			     sgn * (intx + rint((x-intx) * p10) / p10)
			 }
			 else { ## digits < 0
			     p10 <- pow10(-digits)
			     sgn * rint(x/p10) * p10
			 }
		     },
		     "signif" = { ## following ~/R/D/r-devel/R/src/nmath/fprec.c :
                         if(missing(digits)) digits <- 6L
			 if(digits > max(.getPrec(x)) * log10(2))
			     return(x)
			 if(digits < 1) digits <- 1L
			 l10 <- log10(x)
			 e10 <- digits - 1L - floor(l10)
			 r <- x
			 pos.e <- (e10 > 0) ##* 10 ^ e, with e >= 1 : exactly representable
			 if(any(pos.e)) {
			     p10 <- pow10(e10[pos.e])
			     r[pos.e] <- sgn[pos.e]* rint(x[pos.e]*p10) / p10
			 }
			 if(any(neg.e <- !pos.e)) {
			     p10 <- pow10(-e10[neg.e])
			     r[neg.e] <- sgn[neg.e]* rint(x[neg.e]/p10) * p10
			 }
			 r
		     },
		     stop(gettextf("Non-Math2 group generic '%s' -- should not happen",
				   .Generic)))
	  })

##---- mpfrArray / mpfrMatrix --- methods -----------------

## not many needed: "mpfrArray" contain "mpfr",
## i.e., if the above methods are written "general enough", they apply directly

setMethod("sign", "mpfrArray",
	  function(x) structure(.mpfr.sign(x),
				dim = dim(x),
				dimnames = dimnames(x)))
