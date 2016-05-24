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


## error message for invalid seed
.invalid.seed.msg <-
    "invalid 'seed' [NULL / integer(n.sim) / list / NA / \"seq\"]"

##' @title Create a List of Advanced .Random.seed's for "L'Ecuyer-CMRG"
##' @param n number of steps to advance .Random.seed
##' @return a list of length n containing the advanced .Random.seed's
##' @author Marius Hofert and Martin Maechler
LEseeds <- function(n) {
    stopifnot(identical(RNGkind()[1], "L'Ecuyer-CMRG"), n >= 1)
    r <- vector("list", n)
    r[[1]] <- .Random.seed
    for(i in seq_len(n-1))
	r[[i+1]] <- nextRNGStream(r[[i]]) # from parallel
    r
}

## printInfo
printInfo <- local({

    strv <- function(val)
	paste("value:",
              sub("^ *", " ",capture.output(str(unname(val)))))

    SysI <- function()
        paste(format(Sys.time(), "%H:%M:%S"), "on", Sys.info()[["nodename"]])

    F <- function(cc) format(cc, justify = "right")

    formG <- function(gr, j) {
	paste(names(gr), as.matrix(F(gr))[j,], sep= "=", collapse= ", ")
    }

    ##' @title Print Diagnostics for subjob()
    ##' @param i.sim integer running in 1:n.sim
    ##' @param j integer giving the row of the physical grid pGrid
    ##'        (i.sim and j together determine the subjob() argument i)
    ##' @param pGrid physical grid; see subjob()
    ##' @param res4 4-list containing value, error, warning, and time
    ##' @param n.sim: "the" n.sim, possibly NULL
    ##' @return string with information about the sub-job just finished
    ##' @author Marius Hofert and Martin Maechler
    defMoni <- function(i.sim, j, pGrid, res4, n.sim, file="") {
	cat(SysI(), if(length(n.sim) && n.sim > 1)
	    sprintf(": i.sim=%*d, ", 1+floor(log10(n.sim)), i.sim) else ": ",
	    formG(pGrid, j), "; ", strv(res4$value), "\n",
	    sep="", file=file, append=(file != ""))
    }

    list(default = defMoni,
         gfile = function(...) defMoni(..., file="subjob-monitor.txt"),
         fileEach = function(i.sim, j, ...) {
             defMoni(i.sim, j, ..., file=sprintf("subjob-mon-%04d.txt", j))
         })
})

##' @title Compute one Row of the Virtual Grid (subjob)
##' @param i row number of the virtual grid
##' @param pGrid physical grid with all combinations of variables of type "grid" as
##'        returned by mkGrid()
##' @param nonGrids values of non-"grid"-variables (if provided, passed to doOne())
##' @param n.sim number of simulation replications
##' @param seed one of:
##'     NULL: .Random.seed remains untouched; if it doesn't exist, generate
##'           it by calling runif(1); non-reproducible
##'     numeric(n.sim): seeds (numbers) for each of the n.sim simulation
##'                     replications (same seed for each row in the (physical)
##'                     grid); reproducible
##'     vector("list", n.sim): seeds (vectors) for each of the n.sim simulation
##'                            replications (same seed for each row in the (physical)
##'                            grid); this case is meant for the rng L'Ecuyer-CMRG;
##'                            reproducible
##'     NA: .Random.seed remains untouched; if it doesn't exist, so be it;
##'         no 5th component is concatenated to the result of the doOne() call;
##'         non-reproducible
##'     character string: specifying a seeding method; currently only "seq" for
##'                       the seeds 1:n.sim for the n.sim simulation replications;
##'                       reproducible
##' @param keepSeed logical indicating whether .Random.seed is 'stored'/appended
##' @param repFirst logical; if TRUE, first all rows of the (physical) grid are computed for a fixed
##'             replication until the next replication is considered.
##' @param doOne function for computing one row in the (physical) grid; must return
##'        the return value of doCallWE()
##' @param monitor logical or function, see help file
##' @param ... additional arguments passed to doOne()
##' @return a vector of length 5 if(keepSeed), otherwise of length 4. ...
##' @note subjob() is an auxiliary function called from doLapply(), doMclapply(),
##'       doClusterApply() etc.
##' @author Marius Hofert and Martin Maechler
##' { subjob }
subjob <- function(i, pGrid, nonGrids, n.sim, seed, keepSeed=FALSE,
                   repFirst=TRUE, doOne,
                   timer=mkTimer(gcFirst=FALSE), monitor=FALSE, ...)
{
    ## i |-> (i.sim, j) :
    ## determine corresponding i.sim and row j in the physical grid
    if(repFirst) {
	i.sim <- 1 + (i-1) %%  n.sim ## == i  when n.sim == 1
	j     <- 1 + (i-1) %/% n.sim ## row of pGrid
	## Note: this case first iterates over i.sim, then over j:
	## (i.sim,j) = (1,1), (2,1), (3,1),..., (1,2), (2,2), (3,2), ...
    } else {
	ngr <- nrow(pGrid) # number of rows of the (physical) grid
	j     <- 1 + (i-1) %%  ngr ## row of pGrid
	i.sim <- 1 + (i-1) %/% ngr
	## Note: this case first iterates over j, then over i.sim:
	## (i.sim,j) = (1,1), (1,2), (1,2),..., (2,1), (2,2), (2,3), ...
    }

    ## seeding
    if(is.null(seed)) {
	if(!exists(".Random.seed")) runif(1) # guarantees that .Random.seed exists
	## => this is typically not reproducible
    }
    else if(is.numeric(seed)) {
	if(length(seed) != n.sim) stop("'seed' has to be of length ", n.sim)
	set.seed(seed[i.sim]) # same seed for all runs within the same i.sim
	## => calculations based on same random numbers as much as possible
    }
    ## else if(length(seed) == n.sim*ngr && is.numeric(seed)) {
    ##    set.seed(seed[i]) # different seed for *every* row of the virtual grid
    ##    always (?) suboptimal (more variance than necessary)
    ## }
    else if(is.list(seed)) { # (currently) L'Ecuyer-CMRG
	if(length(seed) != n.sim) stop("'seed' has to be of length ", n.sim)
	if(!exists(".Random.seed"))
	    stop(".Random.seed does not exist - in l'Ecuyer setting")
	assign(".Random.seed", seed[[i.sim]], envir = globalenv())
    }
    else if(is.na(seed)) {
	keepSeed <- FALSE
    }
    else {
	if(!is.character(seed)) stop(.invalid.seed.msg)
	switch(match.arg(seed, choices = c("seq")),
	       "seq" = { # sequential seed :
		   set.seed(i.sim) #same seed for all runs within the same i.sim
		   ## => calculations based on the same random numbers
	       },
	       stop("invalid character 'seed': ", seed)
	       )
    }
    ## save seed, compute and return result for one row of the virtual grid
    if(keepSeed) rs <- .Random.seed # <- save here in case it is advanced in doOne

    ## monitor checks happen already in caller!
    if(isTRUE(monitor)) monitor <- printInfo[["default"]]

    ## doOne()'s arguments, grids, non-grids, and '...':
    args <- c(pGrid[j, , drop=FALSE],
	      ## [nonGrids is never missing when called from doLapply() etc.]
	      if(missing(nonGrids) || length(nonGrids) == 0)
	      list(...) else c(nonGrids, ...))
    nmOne <- names(formals(doOne))
    if(!identical(nmOne, "...")) {
	io <- match(nmOne, names(args))         ## __ keep in sync with doCheck() below __
	## adjust *order* of arguments 'args' for doOne()
	## MM: why?? is only needed if some args are not named; do we really want to support that?
	## TODO ?? If there is a "..." among nmOne, match all the others ...
	if(anyNA(io)) {
	    if(!identical("...", (na.One <- nmOne[is.na(io)])) &&
               i.sim == 1 && j == 1) { # <- message *once* only
		## No warning: e.g., '...' argument should not warn
		message(paste(
		"FYI: Argument names from doOne() not present in varlist or '...' of subjob().",
                    "\nNon-matching argument names (may be perfectly valid):",
		    paste(sQuote(na.One), collapse=", ")))
	    }
	    io <- io[!is.na(io)]
	}
	if(length(io))
	    args <- args[c(io, seq_along(args)[-io])]
    }

    r4 <- doCallWE(doOne, args, timer = timer)

    ## monitor (after computation)
    if(is.function(monitor)) monitor(i.sim, j=j, pGrid=pGrid, n.sim=n.sim, res4=r4)

    c(r4, if(keepSeed) list(.Random.seed = rs)) # 5th component .Random.seed
}
##' { end }

##' Check a user's "doOne" function; do similar things as mkAL() and/or subjob()
doCheck <- function(doOne, vList, nChks = ng, verbose=TRUE)
{
    stopifnot(is(vList, "varlist"))
    pGrid <- mkGrid(vList)
    nonGrids <- get.nonGrids(vList)$nonGrids
    stopifnot(is.function(doOne), is.data.frame(pGrid), is.list(nonGrids),
	      is.integer(ng <- nrow(pGrid)), ng >= 1)
    j.s <- sample.int(ng, nChks)
    for(j in j.s) {
	if(verbose) cat(sprintf("j =%3d:", j))
	## doOne()'s arguments, grids, non-grids (but no '...' here [change?])
	gr <- pGrid[j, , drop=FALSE]
	args <- if(length(nonGrids) == 0) as.list(gr) else c(gr, nonGrids)
        ## __ keep in sync with subjob() above __
	nmOne <- names(formals(doOne))
	io <- match(nmOne, names(args))
	if(anyNA(io)) {
	    if(!identical("...", (na.One <- nmOne[is.na(io)])) &&
               j == 1) { # <- message *once* only
		message(paste(
		"FYI: Argument names from doOne() not present in varlist.",
                    "\nNon-matching argument names (may be perfectly valid):",
		    paste(sQuote(na.One), collapse=", ")))
	    }
	    io <- io[!is.na(io)]
	}
	if(length(io))
	    args <- args[c(io, seq_along(args)[-io])]

        r <- doCallWE(doOne, args)
        if(is.null(r$error))
            r <- r$value
	else {
	    message("sample ", j," gave error (that we caught): ",
		    r$error$message)
	    next # j in for()
	}
	dr <- dim(r)
	if(verbose) {
	    cat(" -> length(r) = ",length(r)," dim(r):"); print(dr)
	    cat("str(dimnames(r)):"); str(dimnames(r))
        }
	stopifnot(is.numeric(r))
        if(is.null(dr)) dr <- length(r)
	## TODO MM: check "inner" --map--> dim & dimnames of r
        ## -------  in the same way as mkAL()
        dnr <- dimnames(r)
        ni <- length(l.i <- sapply(getEl(vList, "inner"), length))
        if(anyDuplicated(dr) && is.null(names(dnr)))
            warning("some dimensions in your doOne() result are equal\n",
                    "and cannot be distinguished by names(dimnames(.))")
        ## TODO MM: check names(dnr) *iff* mkAL() uses the same, i.e. not yet
        if(!all(tail(dr, ni) == l.i))
	    stop(gettextf("tail(dim(doOne(..)), %d) differs from %s",
			  ni, 'sapply(getEl(vList, "inner"), length)'),
                 domain=NA)
    }
    invisible(r)
}
