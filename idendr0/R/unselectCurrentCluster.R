unselectCurrentCluster<-function
### Unselect the current cluster.
##keyword<<internal
(
    df, ##<< shared data frame
    dbg=FALSE ##<< debug flag/level
) {
    if (dbg) cat('unselectCurrentCluster called\n')
    if (dbg) printVar(df$currentCluster)

    if (!is.null(df$clusters) && length(df$clusters)>=df$currentCluster
        && !is.null(df$clusters[[df$currentCluster]]$indices)) {
        # remember the current selection
        df<-pushSelectionHistory(df,dbg)
        # unselect the current cluster
        currentClusterIdxInH<-max(df$clusters[[df$currentCluster]]$indices)
        currentClusterLeafs<-computeMemberIndices(df$h,currentClusterIdxInH)
        df$leafColorIdxs[currentClusterLeafs]<-0
        df$unselectedBranches$indices<-
            c(df$unselectedBranches$indices,df$clusters[[df$currentCluster]]$indices)
        df$unselectedBranches$branches<-
            rbind(df$unselectedBranches$branches,df$clusters[[df$currentCluster]]$branches)
        df$clusters[[df$currentCluster]]$indices<-NULL
        df$clusters[[df$currentCluster]]$branches<-NULL
        selectionChanged<-TRUE
    } else {
        # the current cluster was selected, nothing to do
        selectionChanged<-FALSE
    }
    return(list(df=df,selectionChanged=selectionChanged))
    ### a list of shared data frame 'df' and a boolean flag
    ### 'selectionChanged' determing if clsuter selection has changed
    ### (so the caller can learn whether to redraw clusters).
}
