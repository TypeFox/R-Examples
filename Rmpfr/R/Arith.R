#### Define mpfr methods for Arith + Compare + Logic	group functions
####			    ======   =======   =====

### "Math" are done in ./Math.R ,  "Summary" in ./Summary.R
###  ----		 ~~~~~~	    -------	  ~~~~~~~~~
### NB:	 Look at /usr/local/app/R/R_local/src/Brobdingnag/R/brob.R
###					      -----------

if(FALSE) {
 print(getGroupMembers("Ops"))#  "Arith"   "Compare" "Logic"
 .Ops.list <- sapply(getGroupMembers("Ops"),
                     getGroupMembers, simplify=FALSE)
 str(.Ops.list, vec.len = 20)
 ## $ Arith  : chr [1:7] "+" "-" "*" "^" "%%" "%/%" "/"
 ## $ Compare: chr [1:6] "==" ">" "<" "!=" "<=" ">="
 ## $ Logic  : chr [1:2] "&" "|"
}

## Using "vector" and "array" seperately, rather than "ANY"
## ===> shorter distance in method dispatch calculation :
setMethod("Ops", signature(e1 = "mpfr", e2 = "vector"),
	  function(e1, e2) callGeneric(e1, as(e2, "numeric")))
setMethod("Ops", signature(e1 = "vector", e2 = "mpfr"),
	  function(e1, e2) callGeneric(as(e1, "numeric"), e2))

## These should not trigger anymore (because we have "Arith"/"Compare"/...):
setMethod("Ops", signature(e1 = "mpfr", e2 = "array"),
	  function(e1, e2)
	      stop(gettextf("'%s'(mpfr,array) method is not implemented yet",
			    .Generic)))
setMethod("Ops", signature(e1 = "array", e2 = "mpfr"),
	  function(e1, e2)
	      stop(gettextf("'%s'(array,mpfr) method is not implemented yet",
			    .Generic)))

setMethod("Ops", signature(e1 = "mpfr", e2 = "bigz"),
	  function(e1, e2) callGeneric(e1, .bigz2mpfr(e2)))
setMethod("Ops", signature(e1 = "bigz", e2 = "mpfr"),
	  function(e1, e2) callGeneric(.bigz2mpfr(e1), e2))

setMethod("Ops", signature(e1 = "mpfr", e2 = "bigq"),
	  function(e1, e2) callGeneric(e1, ..bigq2mpfr(e2, pmax(.getPrec(e1), 128L))))
setMethod("Ops", signature(e1 = "bigq", e2 = "mpfr"),
	  function(e1, e2) callGeneric(..bigq2mpfr(e1, pmax(.getPrec(e1), 128L)), e2))


setMethod("Logic", signature(e1 = "mpfr", e2 = "mpfr"),
	  function(e1, e2) callGeneric(as(e1, "numeric"),
				       as(e2, "numeric")))
setMethod("Logic", signature(e1 = "mpfr", e2 = "numeric"),
	  function(e1, e2) callGeneric(as(e1, "numeric"), e2))
setMethod("Logic", signature(e1 = "numeric", e2 = "mpfr"),
	  function(e1, e2) callGeneric(e1, as(e2, "numeric")))
## FIXME?: probably also need	 <array, mpfrArray> etc


###-- 2) ----------- Arith --------------------------------------------------

## R version, no longer used:
.mpfr.negative.R <- function(x) {
    xD <- getDataPart(x)# << currently [2011] *faster* than  x@Data
    for(i in seq_along(x))
	slot(xD[[i]], "sign", check=FALSE) <- - xD[[i]]@sign
    setDataPart(x, xD, check=FALSE) ## faster than  x@Data <- xD
}
.mpfr.negative <- function(x) .Call(Rmpfr_minus, x)

setMethod("Arith", signature(e1 = "mpfr", e2="missing"),
	  function(e1,e2) {
	    switch(.Generic,
		   "+" = e1,
		   "-" = .mpfr.negative(e1),
		   stop(paste("Unary operator", .Generic,
			      "not defined for \"mpfr\" numbers"))
		   )
	  } )


.Arith.codes <-
    c("+" = 1, "-" = 2, "*" = 3, "^" = 4, "%%" = 5, "%/%" =6, "/" = 7)
storage.mode(.Arith.codes) <- "integer"

setMethod("Arith", signature(e1 = "mpfr", e2 = "mpfr"),
	  function(e1, e2) {
	      new("mpfr", .Call(Arith_mpfr, e1, e2, .Arith.codes[.Generic]))
	  })

setMethod("Arith", signature(e1 = "mpfr", e2 = "integer"),
	  function(e1, e2) {
	      new("mpfr", .Call(Arith_mpfr_i, e1, e2, .Arith.codes[.Generic]))
	  })
setMethod("Arith", signature(e1 = "integer", e2 = "mpfr"),
	  function(e1, e2) {
	      new("mpfr", .Call(Arith_i_mpfr, e1, e2, .Arith.codes[.Generic]))
	  })

setMethod("Arith", signature(e1 = "mpfr", e2 = "numeric"),# not "integer"
	  function(e1, e2) {
	      new("mpfr", .Call(Arith_mpfr_d, e1, e2, .Arith.codes[.Generic]))
	  })
setMethod("Arith", signature(e1 = "numeric", e2 = "mpfr"),# not "integer
	  function(e1, e2) {
	      new("mpfr", .Call(Arith_d_mpfr, e1, e2, .Arith.codes[.Generic]))
	  })


###-- 3) ----------- Compare --------------------------------------------------

.Compare.codes <- c("==" = 1, ">" = 2, "<" = 3, "!=" = 4, "<=" = 5, ">=" =6)
storage.mode(.Compare.codes) <- "integer"
## Define "Reverse" codes such that, e.g.,
##  .Compare.codes[ .Compare.codesRev[">="] ]  |--> "<="
.Compare.codesRev   <- .Compare.codes # names() in same order; indices swapped:
.Compare.codesRev[] <- .Compare.codes[c(1, 3:2, 4, 6:5)]

setMethod("Compare", signature(e1 = "mpfr", e2 = "mpfr"),
	  function(e1, e2) {
	      .Call(Compare_mpfr, e1, e2, .Compare.codes[.Generic])
	  })

setMethod("Compare", signature(e1 = "mpfr", e2 = "integer"),
	  function(e1, e2) {
	      .Call(Compare_mpfr_i, e1, e2, .Compare.codes[.Generic])
	  })

setMethod("Compare", signature(e1 = "mpfr", e2 =  "numeric"),# not "integer"
	  function(e1, e2) {
	      .Call(Compare_mpfr_d, e1, e2, .Compare.codes[.Generic])
	  })

setMethod("Compare", signature(e1 = "integer", e2 = "mpfr"),
	  function(e1, e2) {
	      .Call(Compare_mpfr_i, e2, e1,
		    .Compare.codesRev[.Generic])
	  })

setMethod("Compare", signature(e1 = "numeric", e2 = "mpfr"),
	  function(e1, e2) {
	      .Call(Compare_mpfr_d, e2, e1,
		    .Compare.codesRev[.Generic])
	  })

### -------------- mpfrArray ------------------------

.dimCheck <- function(a, b) {
    da <- dim(a)
    db <- dim(b)
    if(length(da) != length(db) || any(da != db))
	stop(gettextf("Matrices must have same dimensions in %s",
		      deparse(sys.call(sys.parent()))),
	     call. = FALSE)
    da
}

setMethod("Arith", signature(e1 = "mpfrArray", e2 = "mpfrArray"),
	  function(e1, e2) {
	      .dimCheck(e1, e2)
	      ## else: result has identical dimension:
	      e1@.Data[] <- .Call(Arith_mpfr, e1, e2, .Arith.codes[.Generic])
	      e1
	  })

setMethod("Arith", signature(e1 = "mpfrArray", e2 = "mpfr"),
	  function(e1, e2) {
	      if(length(e1) %% length(e2) != 0)
		  stop("length of first argument (array) is not multiple of the second argument's one")
	      ## else: result has dimension from array:
	      e1@.Data[] <- .Call(Arith_mpfr, e1, e2, .Arith.codes[.Generic])
	      e1
	  })


## "macro-like	encapsulation" -- using .Call(<registered>, *) for checks
.Arith.num.mpfr <- function(x,y, FUN) {
    if(is.integer(x))
	.Call(Arith_i_mpfr, x,y, .Arith.codes[FUN])
    else
	.Call(Arith_d_mpfr, x,y, .Arith.codes[FUN])
}

.Arith.mpfr.num <- function(x,y, FUN) {
    if(is.integer(y))
	.Call(Arith_mpfr_i, x,y, .Arith.codes[FUN])
    else
	.Call(Arith_mpfr_d, x,y, .Arith.codes[FUN])
}

.Compare.num.mpfr <- function(x,y, FUN) {
    if(is.integer(x))
	.Call(Compare_mpfr_i, y,x, .Compare.codesRev[FUN])
    else
	.Call(Compare_mpfr_d, y,x, .Compare.codesRev[FUN])
}

.Compare.mpfr.num <- function(x,y, FUN) {
    if(is.integer(y))
	.Call(Compare_mpfr_i, x,y, .Compare.codes[FUN])
    else
	.Call(Compare_mpfr_d, x,y, .Compare.codes[FUN])
}

setMethod("Arith", signature(e1 = "array", e2 = "mpfr"),# incl "mpfrArray"
	  function(e1, e2) {
	      if(e2Arr <- !is.null(dim(e2)))
		  .dimCheck(e1, e2)
	      else if(length(e1) %% length(e2) != 0)
		  stop("length of first argument (array) is not multiple of the second argument's one")

	      if(e2Arr) {
		  e2@.Data[] <- .Arith.num.mpfr(e1, e2, .Generic)
		  e2
	      } else {
		  r <- new("mpfrArray")
		  r@Dim <- dim(e1)
		  if(!is.null(dn <- dimnames(e1)))
		      r@Dimnames <- dn
		  r@.Data <-  .Arith.num.mpfr(e1, e2, .Generic)
		  r
	      }
	  })
setMethod("Arith", signature(e1 = "mpfr", e2 = "array"),# "mpfr" incl "mpfrArray"
	  function(e1, e2) {
	      if(e1Arr <- !is.null(dim(e1)))
		  .dimCheck(e1, e2)
	      else if(length(e2) %% length(e1) != 0)
		  stop("length of second argument (array) is not multiple of the first argument's one")

	      if(e1Arr) {
		  e1@.Data[] <- .Arith.mpfr.num(e1, e2, .Generic)
		  e1
	      } else {
		  r <- new("mpfrArray")
		  r@Dim <- dim(e2)
		  if(!is.null(dn <- dimnames(e2)))
		      r@Dimnames <- dn
		  r@.Data <- .Arith.mpfr.num(e1, e2, .Generic)
		  r
	      }
	  })

setMethod("Arith", signature(e1 = "mpfrArray", e2 = "numeric"),
	  function(e1, e2) {
	      if(length(e1) %% length(e2) != 0)
		  stop("length of first argument (array) is not multiple of the second argument's one")
	      e1@.Data[] <- .Arith.mpfr.num(e1, e2, .Generic)
	      e1
	  })

setMethod("Arith", signature(e1 = "numeric", e2 = "mpfrArray"),
	  function(e1, e2) {
	      if(length(e2) %% length(e1) != 0)
		  stop("length of second argument (array) is not multiple of the first argument's one")
	      e2@.Data[] <- .Arith.num.mpfr(e1, e2, .Generic)
	      e2
	  })

setMethod("Arith", signature(e1 = "mpfr", e2 = "mpfrArray"),
	  function(e1, e2) {
	      if(length(e2) %% length(e1) != 0)
		  stop("length of second argument (array) is not multiple of the first argument's one")
	      e2@.Data[] <- .Call(Arith_mpfr, e1, e2, .Arith.codes[.Generic])
	      e2
	  })



setMethod("Compare", signature(e1 = "mpfrArray", e2 = "mpfr"),
	  function(e1, e2) {
	      if(is.null(dim(e2))) {
		  if(length(e1) %% length(e2) != 0)
		      stop("length of first argument (array) is not multiple of the second argument's one")
	      } else .dimCheck(e1, e2)

	      structure(.Call(Compare_mpfr, e1, e2, .Compare.codes[.Generic]),
			dim = dim(e1),
			dimnames = dimnames(e1))
	  })

setMethod("Compare", signature(e1 = "mpfr", e2 = "mpfrArray"),
	  function(e1, e2) {
	      if(is.null(dim(e1))) {
		  if(length(e2) %% length(e1) != 0)
		      stop("length of second argument (array) is not multiple of the first argument's one")
	      } else .dimCheck(e1, e2)

	      structure(.Call(Compare_mpfr, e1, e2, .Compare.codes[.Generic]),
			dim = dim(e2),
			dimnames = dimnames(e2))
	  })

setMethod("Compare", signature(e1 = "mpfr", e2 = "array"),# "mpfr" incl "mpfrArray"
	  function(e1, e2) {
	      if(is.null(dim(e1))) {
		  if(length(e2) %% length(e1) != 0)
		      stop("length of second argument (array) is not multiple of the first argument's one")
	      } else .dimCheck(e1, e2)

	      structure(.Compare.mpfr.num(e1, e2, .Generic),
			dim = dim(e2),
			dimnames = dimnames(e2))
	  })

setMethod("Compare", signature(e1 = "array", e2 = "mpfr"),# incl "mpfrArray"
	  function(e1, e2) {
	      if(is.null(dim(e2))) {
		  if(length(e1) %% length(e2) != 0)
		      stop("length of first argument (array) is not multiple of the second argument's one")
	      } else .dimCheck(e1, e2)

	      structure(.Compare.num.mpfr(e1, e2, .Generic),
			dim = dim(e1),
			dimnames = dimnames(e1))
	  })

setMethod("Compare", signature(e1 = "mpfrArray", e2 = "numeric"),# incl integer
	  function(e1, e2) {
	      if(length(e1) %% length(e2) != 0)
		  stop("length of first argument (array) is not multiple of the second argument's one")
	      structure(.Compare.mpfr.num(e1, e2, .Generic),
			dim = dim(e1),
			dimnames = dimnames(e1))
	  })


setMethod("Compare", signature(e1 = "numeric", e2 = "mpfrArray"),
	  function(e1, e2) {

	      if(length(e2) %% length(e1) != 0)
		  stop("length of second argument (array) is not multiple of the first argument's one")
	      structure(.Compare.num.mpfr(e1, e2, .Generic),
			dim = dim(e2),
			dimnames = dimnames(e2))
	  })
