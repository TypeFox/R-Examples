smoothHeatmapAccordingToClusters<-function
### Smooth heatmap by averaging the observations associated with the
### currently selected clusters. Observations that are not contained
### in any cluster currently selected stay untouched.
##keyword<<internal
(
    df, ##<< shared data frame

    dbg.heatmap.smooth = 0 ##<< debug verbosity level
) {
    if (dbg.heatmap.smooth) cat('recomputing smoothed heatmap\n')
    x<-df$x
    for (i in 1:length(df$clusters)) {
        if (!is.null(df$clusters[[i]]) && length(df$clusters[[i]]$indices)>0) {
            if (dbg.heatmap.smooth) cat(sprintf('  cluster idx %d',i))
            clusterIdxInH<-max(df$clusters[[i]]$indices)
            if (dbg.heatmap.smooth) cat(sprintf('  -> cluster %d\n',clusterIdxInH))
            clusterMembers<-computeMemberIndices(df$h,clusterIdxInH)
            x[clusterMembers,]<-rep(colMeans(x[clusterMembers,,drop=FALSE],na.rm=T),each=length(clusterMembers))
        }
    }
    df$xOrderedSmoothed<-x[df$leafOrder,,drop=F]
    df
}
