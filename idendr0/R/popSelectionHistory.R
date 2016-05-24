popSelectionHistory <- function
### Restore (and discard from history) the last cluster selection.
##keyword<<internal
(
    df, ##<< shared data frame
    dbg=FALSE ##<< debug flag/level
) {
    if (dbg) cat('popSelectionHistory called\n')

    if (length(df$selectionHistory)>0) {
        selection<-df$selectionHistory[[1]]
        df$selectionHistory<-df$selectionHistory[-1]
        df$clusters<-selection$clusters
        df$leafColorIdxs<-selection$leafColorIdxs
        df$unselectedBranches<-selection$unselectedBranches
    } else {
        selection<-NULL
    }

    if (dbg) printVar(length(selection$clusters))
    if (dbg>1) printVar(selection)
    if (dbg) printVar(length(df$selectionHistory))

    if (is.null(selection)) {
        rv<-NULL
    } else {
        rv<-df
    }
    return(rv)
    ### shared data frame holding cluster selection (or NULL if no
    ### selection found in the selection history)
}
