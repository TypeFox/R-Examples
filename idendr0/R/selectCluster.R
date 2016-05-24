selectCluster<-function
### Select the nearest cluster to the position 'pos'.
##keyword<<internal
(
    pos, ##<< (x,y) position in the dendro figure
    df, ##<< shared data frame
    dendroZoom, ##<< the current dendro zoom region (used to determine
    ## the scaling)
    dbg=FALSE ##<< debug flag/level
) {
    if (dbg) print('selectCluster called')

    # remember the current selection
    dfOrig<-df
    # and push it
    df<-pushSelectionHistory(df,dbg)

    if (dbg) print(pos)
    x<-as.numeric(pos)[1]
    y<-as.numeric(pos)[2]
    if (dbg) printVar(x)
    if (dbg) printVar(y)
    gw<-xy2gw(list(x,y))
    if (dbg) printVar(gw)

    ## find nearest cluster merging points
    zoom.limits.gw<-dendro2fig(dendroZoom)
    if (dbg) printVar(zoom.limits.gw)
    # branching points
    if (dbg) printVar(dendro2fig(list(df$h$height,NULL)))
    dendro.g<-df$h$height
    if (df$doFlipG) dendro.g<-df$h$height[df$clusterCount]-dendro.g
    branches.g<-dendro2fig(list(dendro.g,NULL))$g
    if (dbg) printVar(branches.g)
    branches.w<-df$branchCenterOffsets
    if (dbg) printVar(branches.w)
    # distances (measured in the current zoom region)
    d<-((branches.g-gw$g)/diff(zoom.limits.gw$g))^2+((branches.w-gw$w)/diff(zoom.limits.gw$w))^2
    if (dbg) printVar(d)
    # index of the nearest cluster
    minI<-which.min(d)
    if (dbg) printVar(minI)

    # subclusters
    subclusters<-computeSubclusterIndices(df$h,minI)
    if (dbg) printVar(subclusters)
    idx<-clusterId2SegmentIds(subclusters)

    # mark selected dendrogram segments
    df$clusters[[df$currentCluster]] <-
        list(indices=subclusters,
            branches=with(df$allBranches$branches,data.frame(x1s=x1s[idx],x2s=x2s[idx],y1s=y1s[idx],y2s=y2s[idx])))

    # unmark deselected subclusters
    unselectedIndices<-setdiff(1:df$clusterCount,subclusters)
    if (dbg) printVar(unselectedIndices)

    # remove the branches selected from the other clusters
    for (i in seq(along=df$clusters)) {
        if (i!=df$currentCluster && !is.null(df$clusters[[i]])) {
            if (length(intersect(subclusters,df$clusters[[i]]$indices))>0) {
                df$clusters[[i]]$indices<-NULL
                df$clusters[[i]]$branches<-NULL
            } else {
                unselectedIndices<-setdiff(unselectedIndices,df$clusters[[i]]$indices)
            }
        }
    }

    # remove branches of selected clusters from unselectedBranches
    if (dbg) printVar(unselectedIndices)
    df$unselectedBranches$indices<-unselectedIndices
    idx<-clusterId2SegmentIds(unselectedIndices)
    df$unselectedBranches$branches<-with(df$allBranches$branches,data.frame(x1s=x1s[idx],x2s=x2s[idx],y1s=y1s[idx],y2s=y2s[idx]))

    # compute leaf colors
    df$leafColorIdxs<-computeLeafColorIdxs(df)

    if (all(df$leafColorIdxs==dfOrig$leafColorIdxs)) {
      # selection has not changed, stick to the previous one
      # (it includes disregarding the push operation)
      df<-dfOrig
    }

    if (dbg>1) printVar(df)
    return(df)
    ### shared data frame
}
