## erf(), erfc()

erf <- function(x) {
    if(is.numeric(x)) 2 * pnorm(x * sqrt(2)) - 1
    else if(is(x, "mpfr")) { # maybe also mpfrMatrix
	##new("mpfr", .Call(Math_mpfr, x, .Math.codes[["erf"]]))
	x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["erf"]])
	x
    }
    else stop("invalid class(x): ", class(x))
}
##    pnorm(x* sqrt(2)) = (1 + erf(x))/2
##==> pnorm(x.)	 = (1 + erf(x./sqrt(2)))/2

##    pnorm(x* sqrt(2), lower=FALSE) = erfc(x)/2
##==> pnorm(x., lower=TRUE)  = erfc(x./sqrt(2))/2
erfc <- function(x) {
    if(is.numeric(x)) 2 * pnorm(x * sqrt(2), lower.tail = FALSE)
    else if(is(x, "mpfr")) {
	x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["erfc"]])
	x
    }
    else stop("invalid class(x): ", class(x))
}

pnorm <- function (q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
{
    if(is.numeric(q) && is.numeric(mean) && is.numeric(sd))
	stats__pnorm(q, mean, sd, lower.tail=lower.tail, log.p=log.p)
    else if((q.mp <- is(q, "mpfr")) || is(mean, "mpfr") || is(sd, "mpfr")) {
	stopifnot(length(lower.tail) == 1, length(log.p) == 1)
	rr <- q <- ((if(q.mp) q else as(q, "mpfr")) - mean) / sd
	if(any(neg <- (q < 0))) ## swap those:	Phi(-z) = 1 - Phi(z)
	    rr[neg] <- pnorm(-q[neg], lower.tail = !lower.tail, log.p=log.p)
	if(any(pos <- !neg)) {
	    q <- q[pos]
	    prec.q <- max(.getPrec(q))
	    rt2 <- sqrt(mpfr(2, prec.q))
	    rr[pos] <- if(lower.tail) {
		eq2 <- erf(q/rt2)
		if(log.p && any(sml <- abs(eq2) < .5)) {
		    r <- q
		    r[ sml] <- log1p(eq2[sml]) - log(2)
		    r[!sml] <- log((1 + eq2[!sml])/2)
		    r
		}
		else {
		    r <- (1 + eq2)/2
		    if(log.p) log(r) else r
		}
	    } else { ## upper.tail
		r <- erfc(q/rt2) / 2
		if(log.p) log(r) else r
	    }
	}
	rr
    } else stop("(q,mean,sd) must be numeric or \"mpfr\"")
}#{pnorm}

dnorm <- function (x, mean = 0, sd = 1, log = FALSE) {
    if(is.numeric(x) && is.numeric(mean) && is.numeric(sd))
	stats__dnorm(x, mean, sd, log=log)
    else if((x.mp <- is(x, "mpfr")) || is(mean, "mpfr") || (s.mp <- is(sd, "mpfr"))) {
	## stopifnot(length(log) == 1)
	prec <- pmax(53, getPrec(x), getPrec(mean), getPrec(sd))
	if(!x.mp) x <- mpfr(x, prec)
	x <- (x - mean) / sd
	twopi <- 2*Const("pi", prec)
	if(!s.mp) sd <- mpfr(sd, prec)
	## f(x) =  1/(sigma*sqrt(2pi)) * exp(-1/2 x^2)
	if(log) ## log( f(x) ) = -[ log(sigma) + log(2pi)/2 + x^2 / 2]
	    -(log(sd) + (log(twopi) + x*x)/2)
	else exp(-x^2/2) / (sd*sqrt(twopi))
    } else stop("invalid arguments (x,mean,sd)")
}

dpois <- function (x, lambda, log = FALSE) {
    if(is.numeric(x) && is.numeric(lambda)) ## standard R
	stats__dpois(x, lambda, log=log)
    else if((l.mp <- is(lambda, "mpfr")) || (x.mp <- is(x, "mpfr"))) {
	prec <- pmax(53, getPrec(lambda), getPrec(x))
	if(!l.mp) lambda <- mpfr(lambda, prec)
	if(!x.mp) x <- mpfr(x, prec)
	if(log)	 -lambda  + x*log(lambda) - lfactorial(x)
	else exp(-lambda) * lambda^x	  /  factorial(x)
    } else
	stop("(x,lambda) must be numeric or \"mpfr\"")
}

dbinom <- function (x, size, prob, log = FALSE) {
    if(is.numeric(x) && is.numeric(size) && is.numeric(prob)) ## standard R
	stats__dbinom(x, size, prob, log=log)
    else if((s.mp <- is(size, "mpfr")) |
	    (p.mp <- is(prob, "mpfr")) |
	    (x.mp <- is(x,    "mpfr"))) {
	stopifnot(is.whole(x)) # R's dbinom() gives NaN's with a warning..
	prec <- pmax(53, getPrec(size), getPrec(prob), getPrec(x))
	if(!s.mp) size <- mpfr(size, prec)
	if(!p.mp) prob <- mpfr(prob, prec)
	if(!x.mp)    x <- mpfr(x, prec)
	## n:= size, p:= prob,	compute	 P(x) = choose(n, x) p^x (1-p)^(n-x)
	C.nx <- chooseMpfr(size, x)
	if(log) log(C.nx) + x*log(prob) + (size-x)*log1p(-prob)
	else C.nx * prob^x * (1-prob)^(size-x)
    } else
	stop("(x,size, prob) must be numeric or \"mpfr\"")
}## {dbinom}



## zeta()
zeta <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr") # keep "mpfrArray"
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["zeta"]])
    x
}

## "FIXME" -- rather use 'bigq' in gmp and the "sumBin" algorithm from copula!
Bernoulli <- function(k, precBits = 128) {
    ## Purpose: Bernoulli Numbers (in high precision)
    ## -----------------------------------------------------------
    ## Arguments: k: non-negative integer vector
    ## -----------------------------------------------------------
    ## Author: Martin Maechler, Date: 12 Dec 2008, 11:35
    stopifnot(all(k >= 0), k == as.integer(k))
    r <- - k * zeta(if(is(k, "mpfr")) 1 - k else mpfr(1 - k, precBits=precBits))
    if(any(k0 <- k == 0)) r[k0] <- mpfr(1, precBits=precBits)
    r
}


## eint() "Exponential integral"
Ei <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr") # keep "mpfrArray"
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["Eint"]])
    x
}

## Li_2() the dilogarithm
Li2 <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr") # keep "mpfrArray"
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["Li2"]])
    x
}


### ------------- Bessel: ---------
## j0, j1, jn
## y0, y1, yn
j0 <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr") # keep "mpfrArray"
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["j0"]])
    x
}
j1 <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["j1"]])
    x
}
y0 <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["y0"]])
    x
}
y1 <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["y1"]])
    x
}

Ai <- function(x) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(Math_mpfr, x, .Math.codes[["Ai"]])
    x
}

jn <- function(n, x, rnd.mode = c('N','D','U','Z','A')) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(R_mpfr_jn, x, n, match.arg(rnd.mode))
    x
}
yn <- function(n, x, rnd.mode = c('N','D','U','Z','A')) {
    if(!inherits(x, "mpfr")) x <- as(x, "mpfr")
    x@.Data[] <- .Call(R_mpfr_yn, x, n, match.arg(rnd.mode))
    x
}


###-------- 2-argument cases -------

## We want to automatically construct the methods needed:
## But atan2() as argument list and  signature	(y, x)
## where  beta() and lbeta()  have  (a,b) --> cannot treat them identically;
## and treat atan2() speparately

## NB: atan2(), beta() and lbeta() all have implicitGeneric()s in methods with no '...'
## ==  ---> canNOT have 3rd argument : rnd.mode = c('N','D','U','Z','A')
##     ---> using   "N"  instead of    match.arg(rnd.mode)
setMethod("atan2", signature(y = "mpfr", x = "mpfr"),
          function(y, x) new("mpfr", .Call(R_mpfr_atan2, y, x, "N")))
setMethod("atan2", signature(y = "mpfr", x = "numeric"),
          function(y, x) new("mpfr", .Call(R_mpfr_atan2, y, .mpfr(x, 128L), "N")))
setMethod("atan2", signature(y = "numeric", x = "mpfr"),
          function(y, x) new("mpfr", .Call(R_mpfr_atan2, .mpfr(y, 128L), x, "N")))
setMethod("atan2", signature(y = "mpfr", x = "ANY"),
          function(y, x) new("mpfr", .Call(R_mpfr_atan2, y, as(x, "mpfr"), "N")))
setMethod("atan2", signature(y = "ANY", x = "mpfr"),
          function(y, x) new("mpfr", .Call(R_mpfr_atan2, as(y, "mpfr"), x, "N")))

setMethod("atan2", signature(y = "mpfrArray", x = "mpfrArray"),
          function(y, x) {
              if(dim(x) != dim(y))
                  stop("array dimensions differ")
              x@.Data[] <- .Call(R_mpfr_atan2, y, x, "N")
              x
          })
setMethod("atan2", signature(y = "mpfrArray", x = "ANY"),
          function(y, x) {
              if(length(y) %% length(x) != 0)
                  stop("length of first argument (array) is not multiple of the second argument's one")
              y@.Data[] <- .Call(R_mpfr_atan2, y, if(is.numeric(x))
                  .mpfr(x, 128L) else as(x, "mpfr"), "N")
              y
          })
setMethod("atan2", signature(y = "ANY", x = "mpfrArray"),
          function(y, x) {
              if(length(x) %% length(y) != 0)
                  stop("length of second argument (array) is not multiple of the first argument's one")
              x@.Data[] <- .Call(R_mpfr_atan2, if(is.numeric(y))
                  .mpfr(y, 128L) else as(y, "mpfr"), x, "N")
              x
          })

## Using  "macro"  {instead of previous aux. function  mpfrMath2setMeth.a.b() :
for (ff in list(c("beta",  "R_mpfr_beta"),
                c("lbeta", "R_mpfr_lbeta"))) eval(substitute(
{
    setMethod(fname, signature(a = "mpfr", b = "mpfr"),
	      function(a, b) new("mpfr", .Call(Csub, a, b, "N")))
    setMethod(fname, signature(a = "mpfr", b = "numeric"),
	      function(a, b) new("mpfr", .Call(Csub, a, .mpfr(b, 128L), "N")))
    setMethod(fname, signature(a = "numeric", b = "mpfr"),
	      function(a, b) new("mpfr", .Call(Csub, .mpfr(a, 128L), b, "N")))
    setMethod(fname, signature(a = "mpfr", b = "ANY"),
	      function(a, b) new("mpfr", .Call(Csub, a, as(b, "mpfr"), "N")))
    setMethod(fname, signature(a = "ANY", b = "mpfr"),
	      function(a, b) new("mpfr", .Call(Csub, as(a, "mpfr"), b, "N")))

    setMethod(fname, signature(a = "mpfrArray", b = "mpfrArray"),
	      function(a, b) {
		  if(dim(b) != dim(a))
		      stop("array dimensions differ")
		  b@.Data[] <- .Call(Csub, a, b, "N")
		  b
	      })
    setMethod(fname, signature(a = "mpfrArray", b = "ANY"),
	      function(a, b) {
		  if(length(a) %% length(b) != 0)
		      stop("length of first argument (array) is not multiple of the second argument's one")
		  a@.Data[] <- .Call(Csub, a, if(is.numeric(b))
				     .mpfr(b, 128L) else as(b, "mpfr"), "N")
		  a
	      })
    setMethod(fname, signature(a = "ANY", b = "mpfrArray"),
	      function(a, b) {
		  if(length(b) %% length(a) != 0)
		      stop("length of second argument (array) is not multiple of the first argument's one")
		  b@.Data[] <- .Call(Csub, if(is.numeric(a))
				     .mpfr(a, 128L) else as(a, "mpfr"), b, "N")
		  b
	      })

}, list(fname = ff[[1]], Csub = as.symbol(ff[[2]]))))


## hypot()
hypot <- function(x,y, rnd.mode = c('N','D','U','Z','A')) {
    if(is(x, "mpfrArray") || is.array(x)) {
	if(is.array(x)) x <- mpfrArray(x, 128L, dim=dim(x), dimnames(x))
	if(is.array(y)) y <- mpfrArray(y, 128L, dim=dim(y), dimnames(y))
	if(is(y, "mpfrArray")) {
	    if(dim(x) != dim(y))
		stop("array dimensions differ")
	    x@.Data[] <- .Call(R_mpfr_hypot, x, y, match.arg(rnd.mode))
	    x
	} else { ## y is not (mpfr)Array
	    if(length(x) %% length(y) != 0)
		stop("length of first argument (array) is not multiple of the second argument's one")
	    x@.Data[] <- .Call(R_mpfr_hypot, x, as(y, "mpfr"), match.arg(rnd.mode))
	    x
	}
    } else if(is(y, "mpfrArray")) {
	if(length(y) %% length(x) != 0)
	    stop("length of second argument (array) is not multiple of the first argument's one")
	y@.Data[] <- .Call(R_mpfr_hypot, as(x, "mpfr"), y, match.arg(rnd.mode))
	y
    }
    else
	new("mpfr", .Call(R_mpfr_hypot, as(x, "mpfr"), as(y, "mpfr"), match.arg(rnd.mode)))
}

## The Beta(a,b)  Cumulative Probabilities are exactly computable for *integer* a,b:
pbetaI <- function(q, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE,
		   precBits = NULL, rnd.mode = c('N','D','U','Z','A'))
{
    stopifnot(length(shape1) == 1, length(shape2) == 1,
	      is.whole(shape1), is.whole(shape2),
	      shape1 >= 0, shape2 >= 0,
	      length(lower.tail) == 1, length(log.p) == 1,
	      0 <= q, q <= 1, ncp == 0,
	      is.null(precBits) ||
	      (is.numeric(precBits) && is.whole(precBits) && precBits >= 2))
    ## Care for too large (a,b) and "integer overflow".
    ## NB:  below have 0:(b - 1) or 0:(a - 1)
    max.ab <- 2^20
    if(is.na(a <- as.integer(shape1)) || (!lower.tail && a > max.ab))
        stop("a = shape1 is too large for 'lower.tail=FALSE' and the current algorithm")
    if(is.na(b <- as.integer(shape2)) || (lower.tail && b > max.ab))
        stop("b = shape2 is too large for 'lower.tail=TRUE' and the current algorithm")
    n <- a + b - 1L
    pr.x <- getPrec(q, bigq. = 256L)
    if(is.null(precBits)) {
        aq <- abs(as.numeric(q))
        mq <- if(any(po <- aq > 0)) min(aq[po]) else 1 # ==> log = 0
        ## -n*log(|x|): such that 1 - |x|^n does not completely cancel
	precBits <- max(128L, pr.x, -as.numeric(n)*log(mq))
    }
    if(pr.x < precBits || !is(q, "mpfr"))
	q <- mpfr(q, precBits=precBits)

    mpfr1 <- list(.Call(const_asMpfr, 1, 16L, "N")) # as prototype for vapply()
    F <- if(log.p) log else identity

    if(lower.tail) {
	## The prob. is	  P[ X <= x ] = \sum_{k=a}^ n    (n \\ k) x^k (1-x)^(n-k)
        ## but we want to sum from 0 {smallest --> largest} as well:
        ##                P[ X <= x ] = \sum_{k=0}^{b-1} (n \\ k) (1-x)^k x^(n-k)
	k <- 0:(b - 1L)
        FUN.x <- function(x) sum(n.choose.k * (1-x)^k * x^(n-k))
    } else { ## upper tail
	## Prob. P[ X > q ] =  1 - P[X <= q ] = \sum_{k=0}^{a-1} (n \\ k) x^k (1-x)^(n-k)
	k <- 0:(a - 1L)
        FUN.x <- function(x) sum(n.choose.k * x^k * (1-x)^(n-k))
    }
    n.choose.k <- chooseZ(n, k)
    roundMpfr(F(
	## "vapply() for "mpfr"
	new("mpfr",
	    vapply(q, FUN.x, mpfr1))),
	      ## reduce the precision, in order to not "claim wrongly":
	      precBits=precBits, match.arg(rnd.mode))
}
