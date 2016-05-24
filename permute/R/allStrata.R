`allStrata` <- function(n, control)
{
    ## seq vector of observation indices
    v <- seq_len(n)
    ## number of groups
    strata <- getStrata(control, which = "plots")
    lev <- length(levels(strata))
    ## compute nperms on number of levels - for this need Within()
    ## and type == typeP
    type <- getType(control, which = "plots")
    newControl <- how(within = Within(type = type))
    nperms <- numPerms(lev, newControl)
    ## result object
    X <- matrix(nrow = nperms, ncol = length(strata))
    ## store the type
    type <- getType(control, which = "plots")
    mirror <- getMirror(control, which = "plots")
    perms <- if(type == "free") {
        allFree(lev)
    } else if(type == "series") {
        allSeries(lev, nperms = nperms, mirror = mirror)
    } else if(type == "grid") {
        nr <- getRow(control, which = "plots")
        nc <- getCol(control, which = "plots")
        constant <- getConstant(control)
        allGrid(lev, nperms = nperms, nr = nr, nc = nc,
                mirror = mirror, constant = constant)
    } else {
        ## if in here, must have both types == "none"
        ## this is here just in case - need to check if this
        ## is possible given calling function...
        return(v)
    }
    sp <- split(v, strata)
    ## build permutations by permuting the split indices (as list)
    ## then undo the original splitting. This respects original indices
    ## of the samples, even where strata ar not contiguous
    for(i in seq_len(nrow(perms))) {
        X[i, ] <- unsplit(sp[perms[i, ]], strata)
    }
    X
}
