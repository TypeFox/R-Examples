## Copyright (C) 2012-2014 Marius Hofert and Martin Maechler
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 2 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


##' @author Marius Hofert
ul <- function(x) unlist(x, recursive=FALSE)

##' @title Converting a List to an Array of Lists
##' @param x list of length prod(dm) where each element is a list of length 4 (or 5)
##'  containing the named elements "value", "error", "warning", "time" (and
##' ".Random.seed"), the first four as returned by doCallWE()
##' @param vList variable specification list
##' @param repFirst logical
##' @param check logical indicating whether checks are carried out
##' @return array of lists of length 4 (or 5)
##' @author Marius Hofert and Martin Maechler
##' @note Arrays are just vectors and thus have to have elements of length 1.
##'	  Therefore, the result of mkAL()  is accessed with [[,..,]],
##'	  so [[,..,]]$value gives the corresponding 'value', for example.
##' { mkAL }
mkAL <- function(x, vList, repFirst, check=TRUE)
{
    grVars <- getEl(vList, "grid", NA)
    n.sim <- get.n.sim(vList)
    ngr <- prod(vapply(lapply(grVars, `[[`, "value"), length, 1L)) # nrow(pGrid)
    lx <- n.sim * ngr
    if(check) {
	stopifnot(is.list(x))
        if(length(x) != lx)
            stop("varlist-defined grid variable dimensions do not match length(x)")
	if(length(x) >= 1) {
	    x1 <- x[[1]]
	    stopifnot(is.list(x1),
		      c("value", "error", "warning", "time") %in% names(x1))
	}
    }
    if(repFirst) ## reorder x
	x <- x[as.vector(matrix(seq_len(lx), ngr, n.sim, byrow=TRUE))]
    iVals <- getEl(vList, "inner")
    xval <- lapply(x, `[[`, "value")
    iLen <- vapply(iVals, length, 1L)
    n.inVals <- prod(iLen)
    if(check) {
	## vector of all "value" lengths
	v.len <- vapply(xval, length, 1L)
	## NB: will be of length zero, when an error occured !!

	##' is N a true multiple of D? includes equality, but we also true vector
	is.T.mult <- function(N, D) N >= D & {q <- N / D; q == as.integer(q) }

	if(!all(eq <- is.T.mult(v.len, n.inVals))) {
	    ## (!all(len.divides <- v.len %% n.inVals == 0)) {
	    not.err <- vapply(lapply(x, `[[`, "error"), is.null, NA)
	    if(!identical(eq, not.err)) {
		msg <- gettextf(
		  "some \"value\" lengths differ from 'n.inVals'=%d without error",
		  n.inVals)
		if(interactive()) {
		    ## warning() instead of stop():
		    ## had *lots* of computing till here --> want to investigate
		    warning(msg, domain=NA, immediate. = TRUE)
		    cat("You can investigate (v.len, xval, etc) now:\n")
		    browser()
		}
		else stop(msg, domain=NA)
	    }
	    if(all(v.len == 0))
		warning(gettextf(
                  "All \"%s\"s are of length zero. The first error message is\n %s",
                  "value", dQuote(conditionMessage(x[[1]][["error"]]))),
                        domain=NA)
	}
    }

    if(length(iVals) > 0 && length(xval) > 0) {
	## ensure that inner variable names are "attached" to x's "value"s :
	if(noArr <- is.null(di <- dim(xval[[1]])))
	    di <- length(xval[[1]])
	rnk <- length(di)# true dim() induced rank
	nI <- length(iLen)# = number of inner Vars;  iLen are their lengths
	for(i in seq_along(xval)) {
	    n. <- length(xi <- xval[[i]])
	    if(n. == 0) # 'if (check)' above has already ensured this is an "error"
		xi <- NA_real_
	    ## else if (n. != n.inVals)
	    ##	   warning(gettext("x[[%d]] is of wrong length (=%d) instead of %d",
	    ##			   i, n., n.inVals), domain=NA)
	    dn.i <- if(noArr) {
		if(nI == 1) list(names(xi)) else rep.int(list(NULL), nI)
	    } else if(is.null(dd <- dimnames(xi))) rep.int(list(NULL), rnk) else dd
	    ## ==>  rnk := length(di) == length(dn.i)
	    if(rnk == nI)# = length(iVals) = length(iLen)  --  simple matching case
		names(dn.i) <- names(iLen)
	    else { # more complicated as doOne() returned a full vector, matrix ...
                if(rnk != length(dn.i)) warning(
                "dim() rank, i.e., length(dim(.)), does not match dimnames() rank")
		if(nI > rnk) # or rather error?
		    warning("nI=length(iVals) larger than length(<dimnames>)")
		else { # nI<rnk==length(di)==length(dn.i) => find matching dim()
		    ## assume inner variables match the *end* of the array
		    j <- seq_len(rnk - nI)
		    j <- which(di[nI+ j] == iLen[j])
		    if(is.null(names(dn.i))) names(dn.i) <- rep.int("", rnk)
		    names(dn.i)[nI+j] <- names(iLen)[j]
		}
	    }
	    x[[i]][["value"]] <- array(xi, dim=if(noArr)iLen else di, dimnames=dn.i)
	}
    }

    gridNms <- mkNms(grVars, addNms=TRUE)
    dmn <- lapply(gridNms, sub, pattern=".*= *", replacement="")
    dm <- vapply(dmn, length, 1L)
    if(n.sim > 1) {
	dm <- c(dm, n.sim=n.sim)
	dmn <- c(dmn, list(n.sim=NULL))
    }
    ## build array
    array(x, dim=dm, dimnames=dmn)
}
##' { end }

##' @title Converting a List to an Array of Lists, Save and Return It
##' @param x result list from mclapply(), clusterApply() etc.
##' @param vList variable specification list
##' @param sfile file name with extension .rds; if sfile=NULL, then nothing is
##'	   saved (and saveSim() is equal to mkAL())
##' @param check logical indicating whether checks are carried out
##' @param doAL logical indicating whether mkAL() should be called
##' @return 'x' or the array returned by mkAL()
##' @author Martin Maechler
##' { saveSim }
saveSim <- function(x, vList, repFirst, sfile, check=TRUE, doAL=TRUE)
{
    if(doAL) {
	a <- tryCatch(mkAL(x, vList, repFirst=repFirst, check=check),
                      error=function(e) e)
	if(inherits(a, "error")) {
	  warning(paste(
            "Relax..: The simulation result 'x' is being saved;",
            "we had an error in 'mkAL(x, *)' ==> returning 'x' (argument, a list).",
            "  you can investigate mkAL(x, ..) yourself. The mkAL() err.message:",
            conditionMessage(a), sep="\n"))
	  a <- x
	}
    } else a <- x
    if(!is.null(sfile))
	saveRDS(a, file=sfile)
    a
}
##' { end }

##' @title Possibly Read Object from an .rds File
##' @param sfile file name with extension .rds
##' @param msg logical indicating whether a message is printed when an object is read
##' @return the object, or NULL (if the file does not exist)
##' @author Martin Maechler
##' { maybeRead }
maybeRead <- function(sfile, msg=TRUE)
{
    if(is.character(sfile) && file.exists(sfile)) {
	if(msg) message("getting object from ", sfile)
	structure(readRDS(sfile), fromFile = TRUE)
    }
}
##' { end }

##' @title Compute Array of Simulation Result Values
##'  --- HIDDEN auxiliary called from getArray()
##' @param x typically, resulting from mkAL(), an array of 5-lists with
##'	   components "value", "error", "warning", "time", and ".Random.seed".
##' @param err.value numeric value which is used in case of an error
##' @param FUN function to be applied before the array is built
##' @return array of values or err.value (in case of an error)
##' @author Marius Hofert and Martin Maechler
valArray <- function(x, err.value = NA, FUN = NULL)
{
    dmn <- dimnames(x)
    dm <- dim(x)
    rr <- lapply(x, `[[`, "value")
    no.err <- sapply(lapply(x, `[[`, "error"), is.null)
    r <- rr[no.err]
    if(length(r) < 1) stop("no non-error values")
    cdim <- unname(lapply(r, dim))
    all.null <- function(x) do.call(all, lapply(x, is.null))
    if(hasdim <- !all.null(cdim)) {
	if(length(cdim <- unique(cdim)) != 1)
	    stop("\"value\" elements of 'x' have different dimensions")
	cdim <- cdim[[1]]
	cNm <- "dimnames"
	cdmn <- lapply(r, dimnames)
    }
    else { ## e.g., when we have no "inner" and return "just" a vector
	if(any(diff(clen <- vapply(r, length, 1L)) != 0))
	    stop("\"value\" elements of 'x' differ in length")
	cdim <- clen[1]
	cNm <- "names"
	cdmn <- lapply(r, names)
    }
    ## The [c]ommon [d]i[mn]ames {includes names(<vector>)}:
    cdmn <- if(!all.null(cdmn)) {
	if(length(cdmn <- unique(cdmn)) != 1)
	    stop(gettextf("\"value\" elements of 'x' have different %s", cNm),
		 domain=NA)
	## if(length(cdim) > 1) cdmn[[1]] else cdmn
        ## Hmm, better if(..) ?
	if(is.list(cc <- cdmn[[1]])) cc else cdmn
    }## else NULL

    hasdim <- length(cdim) > 1 || cdim > 1
    if((dl <- length(cdim) - length(cdmn)) > 0) # e.g. for unnamed *vector* value
	## add artificial dimnames
	cdmn <- setNames(c(rep.int(list(NULL), dl), cdmn),
			 paste("D", seq_len(dl), sep="."))
    if(is.null(FUN)) {
	FUN <- ul
	if(hasdim) {
	    dm	<- c(cdim, dm)
	    dmn <- c(cdmn, dmn)
	}
    } else stopifnot(is.function(FUN))

    ## stopifnot(  sum(no.err) == length(r) )
    NA.proto <- r[[1]]
    NA.proto[] <- err.value
    rr[!no.err] <- list(NA.proto)
    array(FUN(rr), dim=dm, dimnames=dmn)
}## {valArray}

##' @title Compute Arrays Containing the Simulation Results
##' @param x array of lists with components "value", "error", "warning",
##'	   and "time" as returned by mkAL()
##' @param comp string specifying the component to pick out
##' @param FUN function to be applied after lapply() picks out 'comp' of x
##'	       and before the array is built.
##' { getArray }
getArray <- function(x, comp = c("value", "error", "warning", "time"),
		     FUN = NULL, err.value = NA)
{
    comp <- match.arg(comp)
    if(comp == "value")
	return(valArray(x, err.value=err.value, FUN=FUN))
    ## else :
    dmn <- dimnames(x)
    dm <- dim(x)
    if(is.null(FUN)) {
	FUN <-
	    switch(comp,
		   error =, warning = function(x) !vapply(x, is.null, NA),
		   time = ul)
    } else stopifnot(is.function(FUN))
    array(FUN(lapply(x, `[[`, comp)), dim=dm, dimnames=dmn)
}
##' { end }


##' { array2df }
array2df <- function(x, responseName="value") {
    rk <- length(d <- dim(x))
    ## n.sim needs to be "fixed up" {if it exists at all}:
    if(getRversion() >= "3.1.0")
	as.data.frame.table(x, responseName=responseName,
			    base = list(as.character(seq_len(prod(d[-rk])))))
    else {
	dd <- as.data.frame.table(x, responseName=responseName)
	if("n.sim" %in% names(dd))
	    dd$n.sim <- gl(d[rk], prod(d[-rk]))
	dd
    }
}
##' { end }

if(getRversion() <= "3.0.1") ##
ftable.matrix <- ftable.array <- function(x, ...) ftable(as.table(x), ...)

