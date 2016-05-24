## make sure given DVH list and constraint list have the same IDs / constraints
harmoConstrDVH <-
function(x, constr, byPat=TRUE) {
    UseMethod("harmoConstrDVH")
}

harmoConstrDVH.DVHLst <-
function(x, constr, byPat=TRUE) {
    ## check if constr is not a list (but maybe a data frame)
    if(!(is.list(constr) && !is.data.frame(constr))) {
        if(byPat) {
            IDs     <- x[[1]]$patID
            structs <- list(vapply(x, function(y) y$structure, character(1)))
            structs <- setNames(structs, IDs)
        } else {
            structs <- x[[1]]$structure
            IDs     <- list(vapply(x, function(y) y$patID,     character(1)))
            IDs     <- setNames(IDs, structs)
        }

        constr <- constrDF2L1(constr, byPat=byPat, expand=TRUE,
                              dvhID=IDs, dvhStruct=structs)
    }

    ## subset DVH list
    ## and make sure it has same order (IDs/structures) as constraint list
    if(sum(names(x) %in% names(constr)) < 1L) {
        if(byPat) {
            stop("Selected patient not found")
        } else {
            stop("Selected structure not found")
        }
    }

    sharedNames  <- intersect(names(x), names(constr))
    xShared      <- x[sharedNames]
    constrShared <- constr[sharedNames]

    attr(xShared, which="byPat") <- byPat
    class(xShared) <- "DVHLst"
    
    return(list(x=xShared, constr=constrShared))
}

## TODO: same strategy for IDs / structs as above
harmoConstrDVH.DVHLstLst <-
function(x, constr, byPat=TRUE) {
    ## check if constr is not a list (but maybe a data frame)
    if(!(is.list(constr) && !is.data.frame(constr))) {
        if(byPat) {
            IDs     <- names(x)
            structs <- lapply(x, names)
        } else {
            IDs     <- lapply(x, names)
            structs <- names(x)
        }

        constr <- constrDF2L(constr, byPat=byPat, expand=TRUE,
                             dvhID=IDs, dvhStruct=structs)
    }

    ## subset DVH list of lists
    ## and make sure it has same order (IDs/structures) as constraint list
    if(sum(names(x) %in% names(constr)) < 1L) {
        stop("No selected patients/structures found")
    }

    sharedNames  <- intersect(names(x), names(constr))
    xShared      <- x[sharedNames]
    constrShared <- constr[sharedNames]

    attr(xShared, which="byPat") <- byPat
    class(xShared) <- "DVHLstLst"
    
    return(list(x=xShared, constr=constrShared))
}
