pushDendroZoomHistory <- function
### Save the current dendroZoom (as stored in '.sharedEnv').
##keyword<<internal
(
    df, ##<< shared data frame
    dendroZoom, ##<< dendro zoom to push
    dbg=FALSE ##<< debug flag/level
) {
    if (dbg) cat('pushDendroZoomHistory called\n')
    if (dbg) printVar(dendroZoom)

    df$dendroZoomHistory<-c(list(dendroZoom),df$dendroZoomHistory)

    if (dbg) printVar(length(df$dendroZoomHistory))
    return(df)
    ### shared data frame with dendro zoom pushed on top of the
    ### dendro zoom stack
}
