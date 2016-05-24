popDendroZoomHistory <- function
### Restore (and discard from history) the last dendro zoom.
##keyword<<internal
(
    df, ##<< shared data frame
    dbg=FALSE ##<< debug flag/level
) {
    if (dbg) cat('popDendroZoomHistory called\n')

    if (length(df$dendroZoomHistory)>0) {
        dendroZoom<-df$dendroZoomHistory[[1]]
        df$dendroZoomHistory<-df$dendroZoomHistory[-1]
    } else {
        dendroZoom<-NULL
    }
    if (dbg) printVar(dendroZoom)
    if (dbg) printVar(length(df$dendroZoomHistory))
    rv<-list(df=df,dendroZoom=dendroZoom)
    return(rv)
    ### dendro zoom popped (or NULL if stack empty)
}
