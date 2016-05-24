pushSelectionHistory <- function
### Save the current cluster selection (as stored in df).
##keyword<<internal
(
    df, ##<< shared data frame
    dbg=FALSE ##<< debug flag/level
) {
    if (dbg) cat('pushSelectionHistory called\n')
    if (dbg) printVar(length(df$clusters))
    if (dbg>1) printVar(df$clusters)

    df$selectionHistory<-c(
        list(list(clusters=df$clusters,
            leafColorIdxs=df$leafColorIdxs,
            unselectedBranches=df$unselectedBranches)),
        df$selectionHistory)

    if (dbg) printVar(length(df$selectionHistory))

    return(df)
    ### shared data frame with cluster selection pushed on top of the
    ### selection stack.
}
