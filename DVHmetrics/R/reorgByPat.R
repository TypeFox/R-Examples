## re-organize DVH list of lists
## byPat=TRUE  -> top level is patient
## byPat=FALSE -> top level is structure
reorgByPat <-
function(x, byPat=TRUE) {
    ## flatten list 1 level
    xFlat <- unlist(x, recursive=FALSE)
    xReorg <- if(byPat) {
        ## collect patients
        patIDs <- unique(vapply(xFlat, function(y) y$patID, character(1)))

        ## top level: patID
        xReorg <- lapply(patIDs, function(y) {
            Filter(function(z) { z$patID == y }, xFlat)
        })

        xReorg <- setNames(xReorg, patIDs)

        ## within each patient (top level)
        ## set names and assign DVHLst class
        lapply(xReorg, function(y) {
            attr(y, which="byPat") <- TRUE
            ## collect structures
            structs  <- unique(vapply(y, function(z) z$structure, character(1)))
            class(y) <- "DVHLst"
            setNames(y, structs)
        })
    } else {
        ## collect structures
        structs <- unique(vapply(xFlat, function(y) y$structure, character(1)))

        ## top level: structure
        xReorg <- lapply(structs, function(y) {
            Filter(function(z) { z$structure == y }, xFlat)
        })

        xReorg <- setNames(xReorg, structs)

        ## within each structure (top level)
        ## set names and assign DVHLst class
        lapply(xReorg, function(y) {
            attr(y, which="byPat") <- FALSE
            ## collect patients
            patIDs   <- unique(vapply(y, function(z) z$patID, character(1)))
            class(y) <- "DVHLst"
            setNames(y, patIDs)
        })
    }

    ## add byPat flag and assign DVHLstLst class
    attr(xReorg, which="byPat") <- byPat
    class(xReorg) <- "DVHLstLst"
    xReorg
}
