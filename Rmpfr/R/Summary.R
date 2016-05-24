#### Define mpfr methods for Summary  group functions
####			     =======

### "Math" are done in ./Math.R , "Ops", "Arith", "Logic", "Compare" in ./Arith.R

.Summary.codes <-
    c("max" = 1, "min" = 2, "range" = 3, "prod" = 4, "sum" = 5,
      "any" = 10, "all" = 11)
storage.mode(.Summary.codes) <- "integer"

setMethod("Summary", "mpfr",
	  function(x, ..., na.rm=FALSE) {
	      iop <- .Summary.codes[.Generic]
	      r <- .Call(Summary_mpfr, if(length(x)) c(x, ...) else x, na.rm, iop)
	      if(iop <= 5)
		  new("mpfr", r)
	      else ## any, all :
		  r
	  })


stats__quantile.default <- stats:::quantile.default

setMethod("quantile", "mpfr", # 'mpfr' numbers do not have 'names' slot ... (etc)
	  function(x, ...) {
	      if((match("names", names(list(...)), nomatch = 0L)) == 0L)
		  stats__quantile.default(x, ..., names=FALSE)
	      else ## ... contains 'names = ..'
		  stats__quantile.default(x, ...)
	  })

setMethod("mean", "mpfr", function(x, trim = 0, na.rm = FALSE, ...) {
    if(trim == 0) ## based on sum() :
	sum(x, na.rm=na.rm, ...) / length(x)
    else {
	## cut'n'paste from  mean.default() :
	if (!is.numeric(trim) || length(trim) != 1L || trim < 0)
	    stop("'trim' must be numeric of length one, in  [0, 1/2]")
	if (na.rm)
	    x <- x[!is.na(x)]
	n <- length(x)
	if (anyNA(x))
	    mpfr(NA)
	else if (trim >= 0.5)
	    quantile(x, probs = 0.5, na.rm = FALSE)
	else {
	    lo <- floor(n * trim) + 1
	    hi <- n + 1 - lo
	    mean(sort(x, partial = unique(c(lo, hi)))[lo:hi], na.rm = FALSE)
	}
    }
})

setMethod("median", "mpfr",
	  function(x, na.rm=FALSE) quantile(x, probs = 0.5, na.rm=na.rm))


## FIXME: can do this considerably faster in C:
setMethod("which.max", "mpfr", function(x) which.max(x == max(x)))
setMethod("which.min", "mpfr", function(x) which.max(x == min(x)))
