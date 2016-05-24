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


## This cannot be too strict, as then, <varlist>[i] will fail too often
##  ==> the '##' lines have been deactivated (2013-07-10)
.valid.varlist <- function(object) {
    n <- length(object)
    ## if(!n) return("varlist length is zero")
    nms <- names(object)
    if(length(nms) != n)
	return("names() of the varlist must equal list length")
    if(!all(nzchar(nms)))
	return("names(<varlist>) must be non-empty strings each")
    if(any(vapply(lapply(object, `[[`, "value"), is.null, NA)))
	return("every variable in <varlist> must have a non-empty \"value\" part")
    if(!is.character(types <- vapply(object, `[[`, "", "type")) ||
       !all(nzchar(types)))
	return("every variable in <varlist> must have a valid (string) \"type\"")
    ## if(!any(types == "grid"))
    ##     return("varlist must have at least one \"grid\" typed variable")
    TRUE
}

setClass("varlist", "namedList", validity = .valid.varlist)

varlist <- function(...) {
    n <- length(ls <- list(...))
    stopifnot(nzchar(nms <- names(ls)), length(nms) == n,
              !vapply(values <- lapply(ls, `[[`, "value"), is.null, NA))
    ## create default 'expr' from name :
    ex <- lapply(ls, `[[`, "expr")
    if(any(N <- vapply(ex, is.null, NA)))
        for(i in which(N))
            ls[[i]] <- c(ls[[i]], expr = as.symbol(nms[i]))
    ## create default 'type' from name: "N" if n.sim, "frozen" otherwise
    types <- lapply(ls, `[[`, "type")
    if(any(N <- vapply(types, is.null, NA)))
        for(i in which(N))
            ls[[i]] <- c(ls[[i]], type = if(nms[i] == "n.sim") "N" else "frozen")

    ## every "inner" and every "grid" should have length(value) > 1
    lenV <- vapply(values, length, 1L)
    igr <- types %in% c("grid",  "inner")
    if(any(i1 <- lenV[igr] <= 1))
        warning('"grid" and "inner" variables should have length(value) > 1,\n',
                "But these do not:\n", paste(names(lenV)[igr][i1],
                                             collapse=", "))
    new("varlist", ls)
}

dimnames2varlist <- function(dmn) {
    stopifnot(is.character(nms <- names(dmn)), nzchar(nms),
              "character" == vapply(dmn, class, ""))
    build1 <- function(nm) {
        n <- length(e <- dmn[[nm]])
        list(value = e, type = if(n == 1) "frozen" else "grid")
    }
    do.call(varlist, sapply(nms, build1, simplify=FALSE))
}

`$<-.varlist` <- function(x, name, value) {
    x[[name]] <- value
    new("varlist", x)
}

## "danger": result must remain *valid* varlist after subsetting:
`[.varlist` <- function(x, i, ...) {
    new("varlist", .vl.as.list(x)[i])
}

.vl.as.list <- function(from) {
    r <- getDataPart(from)
    names(r) <- names(from)
    r
}
setAs("varlist","list", .vl.as.list)

print.varlist <- function(x, ...) {
    cat("'varlist' object (extending \"namedList\"), with components\n")
    nx <- names(x)
    fn <- setNames(format(nx), nx)
    Form <- function(EE)
	paste(sub("^list\\(", "(", deparse(EE, control=NULL)), collapse="\n  ")
    for(n in nx)
        cat(fn[[n]], ":", Form(x[[n]]), "\n")
    invisible(x)
}
setMethod(show, "varlist", function(object) print.varlist(object))

##' @title Get Elements of a Variable Specification List
##' @param vList list of variable specifications
##' @param type the variable 'type' to restrict the selection to
##' @param comp either a string containing the component name to pick out or NA
##'             (in which case all components are picked out)
##' @return a named list with the selected components of the variables that match
##'         the provided type
##' @author Martin Maechler
##' @note Get all [El]ements (default: "value") of a specific type from a vList
getEl <- function(vList, type="ALL", comp = "value")
{
    if(length(type) == 0) return(list())
    vl <- .vl.as.list(vList) # == as(vList, "list")
    r <- if(type[[1L]] == "ALL")
        vl else vl[vapply(vl, `[[`, "", "type") %in% type]
    if(is.na(comp)) r else lapply(r, `[[`, comp)
}

##' @title Build Grid from Variable Specification List
##' @param vList list of variable specifications
##' @return a data frame
##' @author Martin Maechler
mkGrid <- function(vList)
{
    stopifnot(is.list(vList), length(vList) == 0 ||
              !is.null(vList[[1]][["value"]]))
    do.call(expand.grid,
	    c(getEl(vList, "grid"),
	      list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)))
}

##' format.default(), also for lists, but then "nicely"
forMat <- function(x, trim = FALSE,
                   justify = c("left", "right", "centre", "none"),
                   addClass = FALSE,
                   ...)
{
    ## if(is.function(x)) x # MH: for functions, do nothing (here)
    ## else
    if(!is.list(x)) format(x, trim=trim, justify=match.arg(justify), ...)
    ## this is better (here) than R 3.0.0's  format.default():
    else if(!is.atomic(x[[1]]) && !is.null(nx <- names(x)) && all(nzchar(nx))) {
	if(addClass)
	    paste(abbreviate(sapply(x, class), minlength=2, dot=TRUE),
		  nx, sep=":")
	else nx
    }
    else { ## is.list
	if (missing(trim))
	    trim <- TRUE
	if (missing(justify))
	    justify <- "none"
	res <- lapply(X = x,
		      function(xx, ...) format.default(
			  ## The next line is different than the default:
			  if(is.list(xx)) unlist(xx) else xx,
			  ...),
		      trim = trim, justify = justify, ...)
	sapply(res, paste, collapse = ", ")
    }
}

##' @title Build Names from Variable Specification List
##' @param vList list of variable specifications
##' @param addNms logical
##' @return a named list of the same \code{\link{length}()} and with
##'         the same \code{\link{names}()} as \code{vList}.
##' @author Martin Maechler
mkNms <- function(vList, addNms=FALSE)
{
    r <- lapply(vList, function(ll)
                forMat(ll[["value"]], trim=TRUE, justify="none"))
    if(addNms) sapply(names(r), function(n)
                      paste0(n,"=", r[[n]]), simplify=FALSE)
    else r
}

##' @title Extract and Possibly Modify n.sim
##' @param vList list of variable specifications
##' @return n.sim if it exists, otherwise 1
##' @author Martin Maechler
get.n.sim <- function(vList) {
    n.sim <- getEl(vList, type="N")$n.sim
    if(is.null(n.sim)) n.sim <- 1 else stopifnot(is.numeric(n.sim), n.sim >= 1)
    n.sim
}

##' @title Set or Change n.sim in a varlist
##' @param vList a varlist
##' @param n a positive integer or NULL
##' @return the modified \code{vList}
##' @author Martin Maechler
set.n.sim <- function(vList, n) {
    if(is.null(n)) { ## remove 'n.sim' from  vList
	vList$n.sim <- NULL
   } else {
	stopifnot(length(n) == 1, n == (n <- round(n)), n >= 1)
	if(is.numeric(vList$n.sim$value))
	    vList$n.sim$value <- n
	else ## "create" an n.sim entry
	    vList$n.sim <- list(value = n, type="N", expr = quote(N[sim]))
    }
    vList
}

##' @title Extract (and Possibly Modify) n.sim and non-"grid"-Variables
##' @param vList list of variable specifications
##' @return list with (modified) n.sim and all variables of type not "grid"
##' @author Martin Maechler
get.nonGrids <- function(vList) {
    type <- unique(vapply(vList, `[[`, character(1), "type"))
    type <- type[type != "grid"]
    nGr <- getEl(vList, type=type)
    n.sim <- nGr$n.sim
    if(is.null(n.sim)) n.sim <- 1 else stopifnot(is.numeric(n.sim), n.sim >= 1)
    nGr$n.sim <- NULL
    list(n.sim=n.sim, nonGrids=nGr)
}

