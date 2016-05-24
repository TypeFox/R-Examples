## Copyright (C) 2012 Marius Hofert and Martin Maechler
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


##' Catching and storing warnings and errors simultaneously
##'
##' Catches and saves both warnings and errors and in the case of
##' a warning, also keeps the computed result.
##'
##' @title Catching and storing warnings and errors simultaneously
##' @param expr assignment or function evaluation
##' @return list with 'value'  : value of expr or simpleError
##'		      'warning': simpleWarning or NULL
##' @author Martin Maechler and Marius Hofert, based on hints from
##' @note Based on Luke Tierney's and Bill Dunlap's suggestions, see
##' \url{https://stat.ethz.ch/pipermail/r-help/2010-December/262626.html}
##' { tryCatch.W.E }
tryCatch.W.E <- function(expr){
    W <- NULL
    w.handler <- function(w){ # warning handler
	W <<- w
	invokeRestart("muffleWarning")
    }
    list(value=withCallingHandlers(tryCatch(expr, error=function(e) e),
	                           warning=w.handler),
	 warning=W)
}
##' { end }

## Not exported, and only used because CRAN checks must be faster and not use > 2 cores.
## use 'export  R_PKG_CHECKING_doExtras=true' to use doExtras for all packages
## or  'export R_simsalapar_check_extra=true' for using doExtras only for simsalapar
## [and unset R_PKG_CHECKING_doExtras etc. to unset]
doXchecks <- function() nzchar(Sys.getenv("R_simsalapar_check_extra")) ||
    identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras"))) # any have to be set by hand (not existing)
## Note: if !doExtras: 'R CMD check simsalapar' and 'R CMD check --as-cran simsalapar'
##       work; if doExtras, only the former works
doExtras <- function() interactive() || doXchecks()

## CRAN request: want to *interactively* debug our package on CRAN machine
## ---- and *never* use > 2 cores ==> cannot use doExtras() for core detection
.parallel.chk.users <- c("hofert", "mhofert", "maechler")
nCores4test <- function() {
    ## Maybe use in the future (R > 3.1.1) ..
    ## dc <- suppressWarnings(detectCores(all.tests = TRUE)) ## all tests -> system() warnings
    dc <- detectCores()
    if(length(dc) != 1 || !is.finite(dc)) dc <- 1L
    if(doXchecks() || ## a "use-all-cores" person who has not set R_PKG_AS_CRAN:
       (is.na(Sys.getenv("R_PKG_AS_CRAN", unset=NA)) && ## 'R_PKG_AS_CRAN' from Martin's script
	Sys.info()[["user"]] %in% .parallel.chk.users))

	dc
    else
	min(2L, dc)
}
