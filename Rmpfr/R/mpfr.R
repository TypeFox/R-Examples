#### All methods for  "mpfr" (and "mpfr1") class
#### apart from coercions and the group methods

setMethod("is.finite", "mpfr",
          function(x) .Call(R_mpfr_is_finite, x))
setMethod("is.infinite", "mpfr",
          function(x) .Call(R_mpfr_is_infinite, x))
## MPFR has only "NaN" ( == "NA"  -- hence these two are identical :
setMethod("is.na", "mpfr",
          function(x) .Call(R_mpfr_is_na, x))
setMethod("is.nan", "mpfr",
          function(x) .Call(R_mpfr_is_na, x))

mpfrIs0 <- function(x) {
    if(is(x, "mpfrArray")) .Call(R_mpfr_is_zero_A, x)
    else .Call(R_mpfr_is_zero, x)
    ## sapply(x, function(.) .@exp == - .Machine$integer.max)
}
mpfr.is.0 <- function(x) {
    .Deprecated("mpfrIs0")
    mpfrIs0(x)
}

.mpfr.is.whole <- function(x) {
    if(is(x, "mpfrArray")) .Call(R_mpfr_is_integer_A, x)
    else .Call(R_mpfr_is_integer, x)
}
mpfr.is.integer <- function(x) {
    .Deprecated(".mpfr.is.whole")
    .mpfr.is.whole(x)
}

## is.whole() is now S3 generic, with default method in gmp
## is.whole <- function(x) {
##     if(is.integer(x) || is.logical(x)) rep.int(TRUE, length(x))
##     else if(is.numeric(x)) x == floor(x)
##     else if(is.complex(x)) x == round(x)
##     else if(is(x,"mpfr")) .mpfr.is.whole(x)
##     else rep.int(FALSE, length(x))
## }
is.whole.mpfr <- function(x) .mpfr.is.whole(x)

## The above for "mpfrArray" :
setMethod("is.finite", "mpfrArray",
	  function(x) .Call(R_mpfr_is_finite_A, x))
setMethod("is.infinite", "mpfrArray",
	  function(x) .Call(R_mpfr_is_infinite_A, x))
## MPFR has only "NaN" ( == "NA"  -- hence these two are identical :
setMethod("is.na", "mpfrArray",
	  function(x) .Call(R_mpfr_is_na_A, x))
setMethod("is.nan", "mpfrArray",
	  function(x) .Call(R_mpfr_is_na_A, x))

mpfr_default_prec <- function(prec) {
    if(missing(prec) || is.null(prec))
	.Call(R_mpfr_get_default_prec)
    else {
	stopifnot((prec <- as.integer(prec[1])) > 0)
	.Call(R_mpfr_set_default_prec, prec)
    }
}

.mpfr.minPrec <- function() .Call(R_mpfr_prec_range, 1L)
.mpfr.maxPrec <- function() .Call(R_mpfr_prec_range, 2L)

## must be sync'ed with enum def. in R_mpfr_get_erange in ../src/utils.c
.erange.codes <- c("Emin", "Emax",
                   "min.emin", "max.emin",
                   "min.emax", "max.emax")
.erange.codes <- setNames(seq_along(.erange.codes), .erange.codes)
## FIXME? better function name ??
.mpfr.erange <- function(kind) {
    stopifnot(length(kind) == 1, is.character(kind))
    if(!any(kind == names(.erange.codes)))
        stop("'kind' must be one of ",
             paste(paste0('"', names(.erange.codes), '"'), collapse=", "))
    .Call(R_mpfr_get_erange, .erange.codes[[kind]])
}


.mpfr.erange.set <- function(kind = c("Emin", "Emax"), value) {
    kind <- match.arg(kind)
    ## value can be double precision, and need be for "64-bit long"
    .Call(R_mpfr_set_erange, .erange.codes[[kind]], value)
}


.mpfr.gmp.numbbits <- function() .Call(R_mpfr_get_GMP_numb_bits)

.mpfrVersion <- function() .Call(R_mpfr_get_version)
mpfrVersion <- function()
    numeric_version(sub("^([0-9]+\\.[0-9]+\\.[0-9]+).*","\\1", .mpfrVersion()))

print.mpfr1 <- function(x, digits = NULL, drop0trailing = TRUE, ...) {
    stopifnot(is(x, "mpfr1"), is.null(digits) || digits >= 1)
    cat("'mpfr1' ",
	format(as(x, "mpfr"), digits=digits, drop0trailing=drop0trailing),
	"\n", sep="")
    invisible(x)
}

setMethod(show, "mpfr1", function(object) print.mpfr1(object))

if(FALSE) ## no longer -- as R CMD check complains about use of non-API R_Outputfile
## For testing, debugging etc
if(.Platform$OS.type != "windows") {## No R_Outputfile (in C) on Windows

.print.mpfr <- function(x, digits = NA, ...) {
    stopifnot(is(x, "mpfr"), is.na(digits) || digits >= 1)
    ## digits = NA --> the inherent precision of x will be used
    if(length(x) >= 1)
	.Call(print_mpfr, x, as.integer(digits))
    invisible(x)
}
}# non-Windows only

## a faster version of getDataPart(.) - as we *KNOW* we have a list
## !! If ever the internal representation of such S4 objects changes, this can break !!
getD <- function(x) { attributes(x) <- NULL; x }

## Get or Set the C-global  'R_mpfr_debug_' variable:
.mpfr.debug <- function(i = NA) .Call(R_mpfr_set_debug, as.integer(i))

print.mpfr <- function(x, digits = NULL, drop0trailing = TRUE,
		       right = TRUE, ...) {
    stopifnot(is(x, "mpfr"), is.null(digits) || digits >= 1)
    ## digits = NA --> the inherent precision of x will be used
    n <- length(x)
    ch.prec <-
	if(n >= 1) {
	    rpr <- range(.getPrec(x))
	    paste("of precision ", rpr[1],
		   if(rpr[1] != rpr[2]) paste("..",rpr[2]), " bits")
	}
    cat(n, "'mpfr'", if(n == 1) "number" else "numbers", ch.prec, "\n")
    if(n >= 1)
	print(format(x, digits=digits, drop0trailing=drop0trailing), ...,
	      right=right, quote = FALSE)
    ## .Call(print_mpfr, x, as.integer(digits))
    invisible(x)
}
setMethod(show, "mpfr", function(object) print.mpfr(object))

## Proposal by Serguei Sokol in order to make  diag() work:
if(FALSE)## << MM is in line with our "as.matrix" methods, but is extreme
setMethod("is.matrix", "mpfr",
          function(x) length(dim(x)) == 2L)
## e.g. M0 <- (M <- cbind(mpfr(1.1, 100)^(98:99)))[,FALSE]; diag(M0)
## gives list() instead of length 0 mpfr

## For matrix indexing:  matrix i |-->  regular i :
.mat2ind <- function(i, dim.x, dimnms.x) {
    ndx <- length(dim.x)
    if(!is.null(di <- dim(i))) {
        if(di[2L] == ndx) {
	    ## k-column Matrix subsetting for array of rank k
	    if(is.character(i)) {
		i <- vapply(seq_along(dim.x), function(j)
			    match(i[,j], dimnms.x[[j]]), seq_len(di[1]))
		if(any(is.na(i)))
		    stop("character matrix index out of limits")
	    }
	    i <- if(is.numeric(i))
		i[,1L] + colSums(t(i[,-1L]-1L)* cumprod(dim.x)[-ndx])
	    else getD(i)
	} else {
	    i <- getD(i)
	}
    }
    i
}

## "[" which also keeps names ... JMC says that names are not support(ed|able)
## ---	for such objects..
.mpfr.subset <- function(x,i,j, ..., drop) {
    nA <- nargs()
    if(nA == 2) { ## x[i] etc -- vector case -- to be fast, need C! --
        ## i <- .mat2ind(i, dim(x), dimnames(x))
        xd <- structure(getD(x)[i], names = names(x)[i])
        if(any(iN <- vapply(xd, is.null, NA))) # e.g. i > length(x)
            xd[iN] <- mpfr(NA, precBits = 2L)
        ## faster than  { x@.Data <- xd ; x }:
        setDataPart(x, xd, check=FALSE)
    } else if(nA == 3 && !is.null(d <- dim(x))) { ## matrix indexing(!)
        ## not keeping dimnames though ...
        message("nargs() == 3	 'mpfr' array indexing ... ")
        new("mpfr", structure(getD(x)[i,j,...,drop=drop], dim = d))
        ## keeping dimnames: maybe try
        ##		     D <- getD(x); dim(D) <- d
        ##		     if(!is.null(dn <- dimnames(x))) dimnames(D) <- dn
        ##		     D <- D[i,,drop=drop]
        ##		     new("mpfr", D)
    }
    else
        stop(gettextf("invalid 'mpfr' subsetting (nargs = %d)",nA))
} ## .mpfr.subset()

.mpfr.msubset <- function(x,i,j, ..., drop) {
    nA <- nargs()
    if(nA == 2) {
        i <- .mat2ind(i, dim(x), dimnames(x))
        xd <- structure(getD(x)[i], names=names(x)[i])
        if(any(iN <- vapply(xd, is.null, NA))) # e.g. i > length(x)
            xd[iN] <- mpfr(NA, precBits = 2L)
        ## faster than  { x@.Data <- xd ; x }:
        setDataPart(x[i], xd, check=FALSE)
    }
    else
        stop(gettext("invalid 'mpfr' matrix subsetting with a matrix (nargs = %d)",nA))
} ## .mpfr.msubset()

### ---------- FIXME:  ./array.R  has other  "mpfrArray" methods for "[" and "[<-" !!!!!!-----------

setMethod("[", signature(x = "mpfr", i = "ANY", j = "missing", drop = "missing"),
          .mpfr.subset)
setMethod("[", signature(x = "mpfrArray", i = "matrix", j = "missing", drop = "missing"),
          .mpfr.msubset)

setMethod("[[", signature(x = "mpfr", i = "ANY"),
	  function(x,i) {
	      if(length(i) > 1L) # give better error message than x@.Data[[i]] would:
		  stop("attempt to select more than one element")
	      xd <- getD(x)[[i]] # also gives error when i is "not ok"
              ## faster than { x@.Data <- list(xd) ; x }
              setDataPart(x, list(xd), check=FALSE)
	  })

## "[<-" :
.mpfr.repl <- function(x, i, ..., value, check = TRUE) {
    if(length(list(...))) ## should no longer happen:
	stop("extra replacement arguments ", deparse(list(...)),
	     " not dealt with")
    ## if(!missing(i)) i <- .mat2ind(i, dim(x), dimnames(x))
    n <- length(xD <- getD(x))
    xD[i] <- value
    if((nn <- length(xD)) > n+1)
	## must "fill" the newly created NULL entries
	xD[setdiff((n+1):(nn-1), i)] <- mpfr(NA, precBits = 2L)
    setDataPart(x, xD, check=check)
}
## FIXME: Should not need this; rather add .mat2ind to .mpfr.repl() above
.mpfr.mrepl <- function(x, i, ..., value, check=TRUE) {
    if(length(list(...))) ## should no longer happen:
	stop("extra replacement arguments ", deparse(list(...)),
	     " not dealt with")
    i <- .mat2ind(i, dim(x), dimnames(x))
    n <- length(xD <- getD(x))
    xD[i] <- value
    if((nn <- length(xD)) > n+1)
	## must "fill" the newly created NULL entries
	xD[setdiff((n+1):(nn-1), i)] <- mpfr(NA, precBits = 2L)
    setDataPart(x, xD, check=check)
}

## value = "mpfr"
setReplaceMethod("[", signature(x = "mpfr", i = "ANY", j = "missing",
				value = "mpfr"),
		 function(x, i, j, ..., value) .mpfr.repl(x, i, ..., value=value))
setReplaceMethod("[", signature(x = "mpfrArray", i = "matrix", j = "missing",
				value = "mpfr"),
		 function(x, i, j, ..., value) .mpfr.mrepl(x, i, ..., value=value))

## for non-"mpfr", i.e. "ANY" 'value', coerce to mpfr with correct prec:
setReplaceMethod("[", signature(x = "mpfr", i = "missing", j = "missing",
				value = "ANY"),
	  function(x,i,j, ..., value)
		 .mpfr.repl(x, , value = mpfr(value, precBits =
				 pmax(getPrec(value), .getPrec(x)))))
setReplaceMethod("[", signature(x = "mpfr", i = "ANY", j = "missing",
				value = "ANY"),
	  function(x,i,j, ..., value) {
	      if(length(xi <- x[i]))
		  .mpfr.repl(x, i, value = mpfr(value, precBits =
				   pmax(getPrec(value), .getPrec(xi))))
	      else x # nothing to replace
	  })
setReplaceMethod("[", signature(x = "mpfrArray", i = "matrix", j = "missing",
				value = "ANY"),
	  function(x,i,j, ..., value) {
	      if(length(xi <- x[i]))
		  .mpfr.mrepl(x, i, value = mpfr(value, precBits =
				    pmax(getPrec(value), .getPrec(xi))))
	      else x # nothing to replace
	  })


## I don't see how I could use setMethod("c", ...)
## but this works "magically"  when the first argument is an mpfr :
c.mpfr <- function(...) new("mpfr", unlist(lapply(list(...), as, Class = "mpfr")))

##  duplicated() now works, checked in ../man/mpfr-class.Rd

## sort() works too  (but could be made faster via faster
## ------  xtfrm() method !  [ TODO ]

setMethod("unique", signature(x = "mpfr", incomparables = "missing"),
	  function(x, incomparables = FALSE, ...)
	  new("mpfr", unique(getD(x), incomparables, ...)))

## This is practically identical to  grid's rep.unit :
rep.mpfr <- function(x, times = 1, length.out = NA, each = 1, ...)
    ## Determine an appropriate index, then call subsetting code
    x[ rep(seq_along(x), times=times, length.out=length.out, each=each) ]



setGeneric("pmin", signature = "...")# -> message about override ...
setGeneric("pmax", signature = "...")

## Check if we should "dispatch" to base
## should be fast, as it should not slow down "base pmin() / pmax()"
## Semantically:  <==> is.atomic(x) && !(is(x, "bigz") || is(x, "bigq"))
pm.ok.base <- function(x, cld = getClassDef(class(x))) is.atomic(x) &&
    (!is.object(x) || { !(extends(cld, "bigz") || extends(cld, "bigq")) })

setMethod("pmin", "mNumber",
	  function(..., na.rm = FALSE) {
	      args <- list(...)
              ## Fast(*) check if "base dispatch" should happen (* "fast" for base cases):
	      ## if((allA <- all(vapply(args, is.atomic, NA))) &&
              ##    ((nonO <- !any(is.obj <- vapply(args, is.object, NA))) ||
              ## {
              ##     cld <- lapply(args, function(.) getClassDef(class(.)))
              ##     cld.o <- cld[is.obj]
              ##     all(vapply(cld.o, extends, NA, "bigz") |
              ##         vapply(cld.o, extends, NA, "bigq")) }))
              if(all(vapply(args, pm.ok.base, NA)))
                  return( base::pmin(..., na.rm = na.rm) )
	      ## else: at least one is "mpfr(Matrix/Array)", "bigz" or "bigq"
	      ## if(!allA || nonO)
              cld <- lapply(args, function(.) getClassDef(class(.)))
              ## else have defined cld above
	      is.m <- vapply(cld, extends, NA, "mpfr")
	      is.q <- vapply(cld, extends, NA, "bigq")
	      is.z <- vapply(cld, extends, NA, "bigz")
	      is.N <- vapply(args, function(x) is.numeric(x) || is.logical(x), NA)
	      if(!any(is.m | is.q | is.z)) # should not be needed -- TODO: "comment out"
		  stop("no \"mpfr\", \"bigz\", or \"bigq\" argument -- wrong method chosen; please report!")
	      N <- max(lengths <- vapply(args, length, 1L))
	      any.m <- any(is.m)
	      any.q <- any(is.q)
	      ## precision needed -- FIXME: should be *vector*
	      mPrec <- max(unlist(lapply(args[is.m], .getPrec)),# not vapply
			   if(any(vapply(args[!is.m], is.double, NA)))
			   .Machine$double.digits,
			   if(any.q) 128L,# arbitrary as in getPrec()
			   unlist(lapply(args[is.z], function(z) frexpZ(z)$exp))# as in getPrec()
			   )
	      ## to be the result :
	      ## r <- mpfr(rep.int(Inf, N), precBits = mPrec)
	      ## more efficient (?): start with the first 'mpfr' argument
	      i.frst.m <- which.max(if(any.m) is.m else if(any.q) is.q else is.z)
	      ## ==> r is "mpfr" if there's any, otherwise "bigq", or "bigz"
	      r <- args[[i.frst.m]]
	      if((n.i <- lengths[i.frst.m]) != N)
		  r <- r[rep(seq_len(n.i), length.out = N)]

	      ## modified from ~/R/D/r-devel/R/src/library/base/R/pmax.R
	      has.na <- FALSE
	      ii <- seq_along(lengths) ## = seq_along(args)
	      ii <- ii[ii != i.frst.m]
	      for(i in ii) {
		  x <- args[[i]]
		  if((n.i <- lengths[i]) != N)
		      x <- x[rep(seq_len(n.i), length.out = N)]
		  n.r <- is.na(r); n.x <- is.na(x)
		  ## mpfr() is relatively expensive
		  if(doM <- any.m && !is.m[i] && !is.N[i]) # "bigz", "bigq"
		      ## r is "mpfr"
		      x <- mpfr(x, precBits = mPrec)
		  else if(doQ <- !any.m && !is.q[i] && !is.N[i]) # "bigz"
		      ## r is "bigq"
		      x <- as.bigq(x)
		  if(has.na || (has.na <- any(n.r, n.x))) {
		      r[n.r] <- x[n.r]
		      x[n.x] <- if(!doM && !doQ) as(r[n.x],class(x)) else r[n.x]
		  }
		  change <- r > x
		  change <- which(change & !is.na(change))
		  r[change] <- x[change]
		  if (has.na && !na.rm)
		      r[n.r | n.x] <- NA
	      }
	      ## wouldn't be ok, e.g for 'bigq' r and args[[1]]:
	      ## mostattributes(r) <- attributes(args[[1L]])
	      ## instead :
	      if(!is.null(d <- dim(args[[1L]]))) dim(r) <- d
	      r
	  })## end { pmin }

setMethod("pmax", "mNumber",
	  function(..., na.rm = FALSE) {
	      args <- list(...)
              ## Fast(*) check if "base dispatch" should happen (* "fast" for base cases):
	      ## if((allA <- all(vapply(args, is.atomic, NA))) &&
              ##    ((nonO <- !any(is.obj <- vapply(args, is.object, NA))) ||
              ## {
              ##     cld <- lapply(args, function(.) getClassDef(class(.)))
              ##     cld.o <- cld[is.obj]
              ##     all(vapply(cld.o, extends, NA, "bigz") |
              ##         vapply(cld.o, extends, NA, "bigq")) }))
              if(all(vapply(args, pm.ok.base, NA)))
                  return( base::pmax(..., na.rm = na.rm) )
	      ## else: at least one is "mpfr(Matrix/Array)", "bigz" or "bigq"
	      ## if(!allA || nonO)
              cld <- lapply(args, function(.) getClassDef(class(.)))
              ## else have defined cld above
	      is.m <- vapply(cld, extends, NA, "mpfr")
	      is.q <- vapply(cld, extends, NA, "bigq")
	      is.z <- vapply(cld, extends, NA, "bigz")
	      is.N <- vapply(args, function(x) is.numeric(x) || is.logical(x), NA)
	      if(!any(is.m | is.q | is.z)) # should not be needed -- TODO: "comment out"
		  stop("no \"mpfr\", \"bigz\", or \"bigq\" argument -- wrong method chosen; please report!")
	      N <- max(lengths <- vapply(args, length, 1L))
	      any.m <- any(is.m)
	      any.q <- any(is.q)
	      ## precision needed -- FIXME: should be *vector*
	      mPrec <- max(unlist(lapply(args[is.m], .getPrec)),# not vapply
			   if(any(vapply(args[!is.m], is.double, NA)))
			   .Machine$double.digits,
			   if(any.q) 128L,# arbitrary as in getPrec()
			   unlist(lapply(args[is.z], function(z) frexpZ(z)$exp))# as in getPrec()
			   )
	      ## to be the result :
	      ## r <- mpfr(rep.int(Inf, N), precBits = mPrec)
	      ## more efficient (?): start with the first 'mpfr' argument
	      i.frst.m <- which.max(if(any.m) is.m else if(any.q) is.q else is.z)
	      ## ==> r is "mpfr" if there's any, otherwise "bigq", or "bigz"
	      r <- args[[i.frst.m]]
	      if((n.i <- lengths[i.frst.m]) != N)
		  r <- r[rep(seq_len(n.i), length.out = N)]

	      ## modified from ~/R/D/r-devel/R/src/library/base/R/pmax.R
	      has.na <- FALSE
	      ii <- seq_along(lengths) ## = seq_along(args)
	      ii <- ii[ii != i.frst.m]
	      for(i in ii) {
		  x <- args[[i]]
		  if((n.i <- lengths[i]) != N)
		      x <- x[rep(seq_len(n.i), length.out = N)]
		  n.r <- is.na(r); n.x <- is.na(x)
		  ## mpfr() is relatively expensive
		  if(doM <- any.m && !is.m[i] && !is.N[i]) # "bigz", "bigq"
		      ## r is "mpfr"
		      x <- mpfr(x, precBits = mPrec)
		  else if(doQ <- !any.m && !is.q[i] && !is.N[i]) # "bigz"
		      ## r is "bigq"
		      x <- as.bigq(x)
		  if(has.na || (has.na <- any(n.r, n.x))) {
		      r[n.r] <- x[n.r]
		      x[n.x] <- if(!doM && !doQ) as(r[n.x],class(x)) else r[n.x]
		  }
		  change <- r < x
		  change <- which(change & !is.na(change))
		  r[change] <- x[change]
		  if (has.na && !na.rm)
		      r[n.r | n.x] <- NA
	      }
	      ## wouldn't be ok, e.g for 'bigq' r and args[[1]]:
	      ## mostattributes(r) <- attributes(args[[1L]])
	      ## instead :
	      if(!is.null(d <- dim(args[[1L]]))) dim(r) <- d
	      r
	  })## end { pmax }


### seq() :

## seq.default()  and  seq.Date()  as examples :
## ~/R/D/r-devel/R/src/library/base/R/seq.R    and
## ~/R/D/r-devel/R/src/library/base/R/dates.R

seqMpfr <- function(from = 1, to = 1, by = ((to - from)/(length.out - 1)),
		    length.out = NULL, along.with = NULL, ...)
{
    if(h.from <- !missing(from)) {
	lf <- length(from)
	if(lf != 1) stop("'from' must be of length 1")
    }
    if (nargs() == 1L && h.from) { # 'One'
	if(is.numeric(from) || is(from,"mpfr")) {
	    to <- from; from <- mpfr(1, getPrec(from))
	} else stop("'from' is neither numeric nor \"mpfr\"")
    }
    ## else if (!is(from, "mpfr")) from <- as(from, "mpfr")

    if(!missing(to)) {
	if (!is(to, "mpfr")) to <- as(to, "mpfr")
	if (length(to) != 1) stop("'to' must be of length 1")
    }
    if (!missing(along.with)) {
	length.out <- length(along.with)
    } else if (!is.null(length.out)) {
	if (length(length.out) != 1) stop("'length.out' must be of length 1")
	length.out <- ceiling(length.out)
    }
##     status <- c(!missing(to), !missing(by), !is.null(length.out))
##     if(sum(status) != 2)
## ## stop("exactly two of 'to', 'by' and 'length.out' / 'along.with' must be specified")
##	   warning("not exactly two of 'to', 'by' and 'length.out' / 'along.with' have been specified")

    if(is.null(length.out)) {
	if(!is(to,   "mpfr")) to   <- as(to,   "mpfr")
	if(!is(from, "mpfr")) from <- as(from, "mpfr")# need it again
	del <- to - from
	if(del == 0 && to == 0) return(to)
	if(missing(by)) {
	    by <- mpfr(sign(del), getD(from)[[1]]@prec)
	}
    }
    if (!is(by, "mpfr")) by <- as(by, "mpfr")
    if (length(by) != 1) stop("'by' must be of length 1")

    ## ---- This is  cut n paste  from seq.default() :
    ## ---- It should work, since "arithmetic works for mpfr :
    if(is.null(length.out)) {
	n <- del/by
	if(!(length(n) && is.finite(n))) {
	    if(length(by) && by == 0 && length(del) && del == 0)
		return(from)
	    stop("invalid (to - from)/by in seq(.)")
	}
	if(n < 0)
	    stop("wrong sign in 'by' argument")
	if(n > .Machine$integer.max)
	    stop("'by' argument is much too small")

	dd <- abs(del)/max(abs(to), abs(from))
	if (dd < 100*.Machine$double.eps) return(from)
	n <- as.integer(n + 1e-7)
	x <- from + (0:n) * by
	## correct for overshot because of fuzz
	if(by > 0) pmin(x, to) else pmax(x, to)
    }
    else if(!is.finite(length.out) || length.out < 0)
	stop("length must be non-negative number")
    else if(length.out == 0)
	as(from,"mpfr")[FALSE] # of same precision
    ## else if (One) 1:length.out
    else if(missing(by)) {
	## if(from == to || length.out < 2) by <- 1
        length.out <- as.integer(length.out)
	if(missing(to))
	    to <- as(from,"mpfr") + length.out - 1
	if(missing(from))
	    from <- to - length.out + 1
	if(length.out > 2)
	    if(from == to)
		rep.int(as(from,"mpfr"), length.out)
	    else { f <- as(from,"mpfr")
		   as.vector(c(f, f + (1:(length.out - 2)) * by, to))
}
	else as.vector(c(as(from,"mpfr"), to))[seq_len(length.out)]
    }
    else if(missing(to))
	as(from,"mpfr") + (0:(as.integer(length.out) - 1L)) * by
    else if(missing(from))
	to - ((as.integer(length.out) - 1L):0) * by
    else stop("too many arguments")
}

if(FALSE) { ##-- --- I don't see *any* way  to define  seq() {S4} methods
    ## 1. Currently  need a  setGeneric() :
    ## ---- just calling setMethod("seq",...) as below fails directly {signature problem}

    ## 2. Trying three different variations --- all of them render the
    ##    *default method invalid :
    ###   --->    seq(1, length.out=3)  # afterwards fails with   " missing 'by' "
setGeneric("seq", function(from, to, by, ...) standardGeneric("seq"),
	   useAsDefault = function(from, to, by, ...)
	   base::seq(from, to, by, ...))

setGeneric("seq", function(from, to, by, ...) standardGeneric("seq"),
	   useAsDefault =
	   function(from = 1, to = 1, by = ((to-from)/(length.out-1)), ...)
	   base::seq(from, to, by, ...))

setGeneric("seq", function (from, to, by, length.out, along.with, ...)
	   standardGeneric("seq"),
	   signature = c("from", "to", "by"),
	   useAsDefault = {
	       function(from = 1, to = 1, by = ((to-from)/(length.out-1)),
			length.out = NULL, along.with = NULL, ...)
		   base::seq(from, to, by,
			     length.out=length.out, along.with=along.with, ...)
	   })

setMethod("seq", c(from = "mpfr", to = "ANY", by = "ANY"), seqMpfr)
setMethod("seq", c(from = "ANY", to = "mpfr", by = "ANY"), seqMpfr)
setMethod("seq", c(from = "ANY", to = "ANY", by = "mpfr"), seqMpfr)

}##--not yet-- defining seq() methods -- as it fails

## the fast mpfr-only version - should *not* return empty, hence the default:
.getPrec <- function(x) {
    if(length(x)) vapply(getD(x), slot, 1L, "prec")
    else mpfr_default_prec()
}

##' The *relevant* number of "bit"/"digit" characters in character vector x
##' (i.e. is vectorized)
.ncharPrec <- function(x, base) {
    if((base ==  2 && any(i <- tolower(substr(x,1L,2L)) == "0b")) ||
       (base == 16 && any(i <- tolower(substr(x,1L,2L)) == "0x"))) {
        i <- which(i)
        x[i] <- substr(x[i], 3L, 1000000L)
    }
    nchar(gsub("[-.]", '', x), "bytes")
}

## the user version
getPrec <- function(x, base = 10, doNumeric = TRUE, is.mpfr = NA, bigq. = 128L) {
    if(isTRUE(is.mpfr) || is(x,"mpfr"))
	vapply(getD(x), slot, 1L, "prec")# possibly of length 0
    else if(is.character(x)) ## number of digits --> number of bits
	ceiling(log2(base) * .ncharPrec(x, base))
    else if(is.logical(x))
	2L # even 1 would suffice - but need 2 (in C ?)
    else if(is.raw(x)) {
	if(is.object(x)) { ## Now deal with 'bigz' and 'bigq'
	    if(inherits(x,"bigz"))
		frexpZ(x)$exp
	    else if(inherits(x,"bigq")) {
		if(missing(bigq.)) {
		    warning("default precision for 'bigq' arbitrarily chosen as ", bigq.)
		    bigq.
		}
		else as.integer(bigq.)
	    }
	    else 8L
	} else 8L
    }
    else {
	if(!doNumeric)
	    stop("must specify 'precBits' for numeric 'x' when 'doNumeric' is false")
	## else
	if(is.integer(x)) 32L
	else if(is.double(x)) 53L
	else if(length(x) == 0) mpfr_default_prec()
	else stop(sprintf("cannot determine 'precBits' for x of type '%s'",
			  typeof(x)))
    }
}


### all.equal()

## TODO ?? <<<<<<<<<<<
## ====
## 2) instead of  as(., "mpfr")	 use  mpfr(., precBits = <smart>)

## For two "mpfr"s, use a  "smart" default tolerance :
setMethod("all.equal", signature(target = "mpfr", current = "mpfr"),
	  function (target, current,
		    tolerance =
		    2^-(0.5 * min(mean(.getPrec(target)),
				  mean(.getPrec(current)))), ...)
      {
	  ## to use "our" mean() :
	  environment(all.equal.numeric) <- environment()
	  all.equal.numeric(target, current, tolerance=tolerance, ...)
      })

setMethod("all.equal", signature(target = "mpfr", current = "ANY"),
	  function (target, current,
		    tolerance = .Machine$double.eps^0.5, ...) {
	      ## to use "our" mean() :
	      environment(all.equal.numeric) <- environment()
	      all.equal.numeric(target, as(current, "mpfr"),
				tolerance=tolerance, ...)
	  })

setMethod("all.equal", signature(target = "ANY", current = "mpfr"),
	  function (target, current,
		    tolerance = .Machine$double.eps^0.5, ...) {
	      ## to use "our" mean() :
	      environment(all.equal.numeric) <- environment()
	      all.equal.numeric(as(target, "mpfr"), current,
				tolerance=tolerance, ...)
	  })

##' This is almost identical to diff.default -- ~/R/D/r-devel/R/src/library/base/R/diff.R
##' But that uses unclass(x) unfortunately
diff.mpfr <- function(x, lag = 1L, differences = 1L, ...)
{
    ismat <- is(x, "mpfrArray") ##_ is.matrix(x)
    xlen <- if(ismat) dim(x)[1L] else length(x)
    if (length(lag) > 1L || length(differences) > 1L ||
        lag < 1L || differences < 1L)
	stop("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen)
	return(x[0L]) # empty, but of proper mode
    i1 <- -seq_len(lag)
    if (ismat)
	for (i in seq_len(differences))
	    x <- x[i1, , drop = FALSE] -
                x[-nrow(x):-(nrow(x)-lag+1L), , drop = FALSE]
    else
        for (i in seq_len(differences))
            x <- x[i1] - x[-length(x):-(length(x)-lag+1L)]
    x
}

str.mpfr <- function(object, nest.lev, give.head = TRUE, digits.d = 12,
                     vec.len = NULL, drop0trailing=TRUE,
                     width = getOption("width"), ...) {
    ## utils:::str.default() gives  "Formal class 'mpfr' [package "Rmpfr"] with 1 slots"
    cl <- class(object)
    le <- length(object)
    if(le == 0) { print(object); return(invisible()) }
    if(isArr <- is(object, "mpfrArray")) di <- dim(object)
    r.pr <- range(getPrec(object))
    onePr <- r.pr[1] == r.pr[2]
    if(give.head)
	cat("Class", " '", paste(cl, collapse = "', '"),
	    "' [package \"", attr(cl, "package"), "\"] of ",
	    if(isArr) paste("dimension", deparse(di, control = NULL))
	    else paste("length", le), " and precision",
	    if(onePr) paste("", r.pr[1]) else paste0("s ", r.pr[1],"..",r.pr[2]),
	    "\n", sep = "")
    if(missing(nest.lev)) nest.lev <- 0
    cat(paste(rep.int(" ", max(0,nest.lev+1)), collapse= ".."))
    ## if object is long, drop the rest which won't be used anyway:
    max.len <- max(100, width %/% 3 + 1, if(is.numeric(vec.len)) vec.len)
    if(le > max.len) object <- object[seq_len(max.len)]
    if(!is.null(digits.d))## reduce digits where precision is smaller:
	digits.d <- pmin(digits.d,
			 ceiling(log(2)/log(10) * .getPrec(object)))
    if(is.null(vec.len)) { # use width and precision (and remain simple enough)
        ff <- formatMpfr(object, digits=digits.d, drop0trailing=drop0trailing, ...)
	nch <- if(getRversion() >= "3.2.1") nchar(ff, keepNA=FALSE) else nchar(ff)
	fits <- !any(too.lrg <- cumsum(nch) + length(nch)-1L > width)
	if(!fits)
	    vec.len <- max(2L, which.max(too.lrg) - 1L)
    } else
	fits <- le <= vec.len
    if(!fits) {
	object <- object[i <- seq_len(vec.len)]
	digits.d <- digits.d[i]
    }
    cat(formatMpfr(object, digits=digits.d, drop0trailing=drop0trailing, ...),
	if(fits) "\n" else "...\n")
} ## {str.mpfr}

