## Extractor functions for blocks, plots and within, plus strata,
## etc ...

## Blocks
`getBlocks` <- function(object, ...) {
    UseMethod("getBlocks")
}

`getBlocks.default` <- function(object, ...) {
    stop("No default method for 'getBlocks()'")
}

`getBlocks.how` <- function(object, ...) {
    object$blocks
}

## Plots
`getPlots` <- function(object, ...) {
    UseMethod("getPlots")
}

`getPlots.default` <- function(object, ...) {
    stop("No default method for 'getPlots()'")
}

`getPlots.how` <- function(object, ...) {
    object$plots
}

## Within plots
`getWithin` <- function(object, ...) {
    UseMethod("getWithin")
}

`getWithin.default` <- function(object, ...) {
    stop("No default method for 'getWithin()'")
}

`getWithin.how` <- function(object, ...) {
    object$within
}

## Strata
`getStrata` <- function(object, ...) {
    UseMethod("getStrata")
}

`getStrata.default` <- function(object, ...) {
    stop("No default method for 'getStrata()'")
}

`getStrata.how` <- function(object,
                                  which = c("plots","blocks"),
                                  drop = TRUE, ...) {
    which <- match.arg(which)
    if(isTRUE(all.equal(which, "plots")))
        strata <- object$plots$strata
    else if(isTRUE(all.equal(which, "blocks")))
        strata <- object$blocks #object$blocks$strata
    else
        stop("Ambiguous `which`")
    if(isTRUE(drop) && !is.null(strata))
        strata <- droplevels(strata)
    strata
}

`getStrata.Plots` <- function(object, drop = TRUE, ... ) {
    strata <- object$strata
    if(isTRUE(drop) && !is.null(strata))
        strata <- droplevels(strata)
    strata
}

## Get type of permutation
`getType` <- function(object, ...) {
    UseMethod("getType")
}

`getType.default` <- function(object, ...) {
    stop("No default method for 'getType()'")
}

`getType.how` <- function(object,
                          which = c("plots","within"), ...) {
    which <- match.arg(which)
  if(isTRUE(all.equal(which, "plots")))
      type <- getPlots(object)$type
  else if(isTRUE(all.equal(which, "within")))
      type <- getWithin(object)$type
  else
      stop("Ambiguous `which`")
  type
}

`getType.Within` <- function(object, ...) {
    object$type
}

`getType.Plots` <- function(object, ...) {
    object$type
}

## suppose we can also have setBlocks() etc...
## to update the control object in place....

## Get mirroring status
`getMirror` <- function(object, ...) {
    UseMethod("getMirror")
}

`getMirror.default` <- function(object, ...) {
    stop("No default method for 'getMirror()'")
}

`getMirror.how` <- function(object,
                                    which = c("plots","within"), ...) {
    which <- match.arg(which)
    if(isTRUE(all.equal(which, "plots")))
        mirror <- getPlots(object)$mirror
    else if(isTRUE(all.equal(which, "within")))
        mirror <- getWithin(object)$mirror
    else
        stop("Ambiguous `which`")
    mirror
}

`getMirror.Within` <- function(object, ...) {
    object$mirror
}

`getMirror.Plots` <- function(object, ...) {
    object$mirror
}

## Get constant status - i.e. same permutation in each Plot
`getConstant` <- function(object, ...) {
    UseMethod("getConstant")
}

`getConstant.default` <- function(object, ...) {
    stop("No default method for 'getConstant()'")
}

`getConstant.how` <- function(object, ...) {
    getWithin(object)$constant
}

`getConstant.Within` <- function(object, ...) {
    object$constant
}

## Get the number of rows and colums from grid designs
`getRow` <- function(object, ...) {
    UseMethod("getRow")
}

`getRow.default` <- function(object, ...) {
    NROW(object)
}

`getRow.how` <- function(object, which = c("plots","within"),
                                 ...) {
    which <- match.arg(which)
    if(isTRUE(all.equal(which, "plots")))
        nrow <- getPlots(object)$nrow
    else if(isTRUE(all.equal(which, "within")))
        nrow <- getWithin(object)$nrow
    else
        stop("Ambiguous `which`")
    nrow
}

`getRow.Within` <- function(object, ...) {
    object$nrow
}

`getRow.Plots` <- function(object, ...) {
    object$nrow
}

`getCol` <- function(object, ...) {
    UseMethod("getCol")
}

`getCol.default` <- function(object, ...) {
    NCOL(object)
}

`getCol.how` <- function(object, which = c("plots","within"),
                                 ...) {
    which <- match.arg(which)
    if(isTRUE(all.equal(which, "plots")))
        ncol <- getPlots(object)$ncol
    else if(isTRUE(all.equal(which, "within")))
        ncol <- getWithin(object)$ncol
    else
        stop("Ambiguous `which`")
    ncol
}

`getCol.Within` <- function(object, ...) {
    object$ncol
}

`getCol.Plots` <- function(object, ...) {
    object$ncol
}

`getDim` <- function(object, ...) {
    UseMethod("getDim")
}

`getDim.default` <- function(object, ...) {
    dim(object)
}

`getDim.how` <- function(object, which = c("plots","within"),
                                 ...) {
    which <- match.arg(which)
    if(isTRUE(all.equal(which, "plots"))) {
        PL <- getPlots(object)
        nc <- PL$ncol
        nr <- PL$nrow
    } else if(isTRUE(all.equal(which, "within"))) {
        WI <- getWithin(object)
        nc <- WI$ncol
        nr <- WI$nrow
    } else {
        stop("Ambiguous `which`")
    }
    c(nr, nc)
}

`getDim.Within` <- function(object, ...) {
    c(object$nrow, object$ncol)
}

`getDim.Plots` <- function(object, ...) {
    c(object$nrow, object$ncol)
}

## return the requested number of permutations
`getNperm` <- function(object, ...) {
    UseMethod("getNperm")
}

`getNperm.default` <- function(object, ...) {
    stop("No default method for `getNperm`")
}

`getNperm.how` <- function(object, ...) {
    object$nperm
}

## Returns maximum permutation threshold
`getMaxperm` <- function(object, ...) {
    UseMethod("getMaxperm")
}

`getMaxperm.default` <- function(object, ...) {
    stop("No default method for `getMaxperm`")
}

`getMaxperm.how` <- function(object, ...) {
    object$maxperm
}

## Returns minimum permutation threshold
`getMinperm` <- function(object, ...) {
    UseMethod("getMinperm")
}

`getMinperm.default` <- function(object, ...) {
    stop("No default method for `getMinperm`")
}

`getMinperm.how` <- function(object, ...) {
    object$minperm
}

## Returns status of complete enumeration
`getComplete` <- function(object, ...) {
    UseMethod("getComplete")
}

`getComplete.default` <- function(object, ...) {
    stop("No default method for `getComplete`")
}

`getComplete.how` <- function(object, ...) {
    object$complete
}

## Returns whether all permutation should/should not be made
`getMake` <- function(object, ...) {
    UseMethod("getMake")
}

`getMake.default` <- function(object, ...) {
    stop("No default method for `getMake`")
}

`getMake.how` <- function(object, ...) {
    object$make
}

## Returns whether the observed permutation should be in
## the set of permutations
`getObserved` <- function(object, ...) {
    UseMethod("getObserved")
}

`getObserved.default` <- function(object, ...) {
    stop("No default method for `getObserved`")
}

`getObserved.how` <- function(object, ...) {
    object$observed
}

## Extractor for all.perms component
`getAllperms` <- function(object, ...) {
    UseMethod("getAllperms")
}

`getAllperms.how` <- function(object, ...) {
    object$all.perms
}

`getAllperms.default` <- function(object, ...) {
    stop("No default method for `getAllperms`")
}
