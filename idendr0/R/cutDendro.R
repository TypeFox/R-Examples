cutDendro<-function
### Select clusters by cutting the current zoom of the dendrogram at a
### specified height.
##keyword<<internal
(
    df, ##<< shared data frame
    cutG, ##<< cutting grow height
    dendroZoom, ##<< current dendro zoom region
    dbg=FALSE ##<< debug level
) {

    if (dbg) printVar(dendroZoom)
    if (dbg) printVar(dendroZoom$g)
    if (dbg) printVar(cutG)
    h<-df$h
    branchCenterOffsets<-df$branchCenterOffsets

    # remember current selection
    df<-pushSelectionHistory(df,dbg)

    if (df$doFlipG) cutG<-df$h$height[df$clusterCount]-cutG
    ch<-cutree(h,h=cutG)
    selectedClusterCandidateCount<-max(ch)
    if (dbg) printVar(selectedClusterCandidateCount)
    if (selectedClusterCandidateCount>1) {
        # allocate new clusters
        newClusters<-vector('list',selectedClusterCandidateCount)
        for (i in 1:selectedClusterCandidateCount) {
            newClusters[[i]]<-list(indices=NULL,branches=NULL)
        }

        parentClusters<-(df$clusterCount-selectedClusterCandidateCount+2):df$clusterCount
        if (dbg>1) printVar(parentClusters)
        selectedClusters<-rep(NA,selectedClusterCandidateCount)
        selectedClustersIdx<-0
        #if (dbg) printVar(currentYlim)
        for (parent in parentClusters) {
            # if left subcluster of the cluster being processed is below
            # the X threshold, take it as a subcluster (otherwise, it is
            # a cluster higher in the hierarchy and its left subcluster is
            # above the X threshold)
            if (!is.element(h$merge[parent,1],parentClusters)) {
                # is the subcluster visible in the current view?
                if (dbg>1) printVar(branchCenterOffsets[h$merge[parent,1]])
                if (branchCenterOffsets[h$merge[parent,1]]>=dendroZoom$w[1] &&
                    branchCenterOffsets[h$merge[parent,1]]<=dendroZoom$w[2]) {
                    if (dbg>1) cat(' selected\n')
                    selectedClustersIdx<-selectedClustersIdx+1
                    selectedClusters[selectedClustersIdx]<-h$merge[parent,1]
                } else if (dbg>1) cat(' not selected\n')
            }
            # the same for the right (2nd) subcluster
            if (!is.element(h$merge[parent,2],parentClusters)) {
                # is the subcluster visible in the current view?
                if (dbg>1) printVar(branchCenterOffsets[h$merge[parent,2]])
                if (branchCenterOffsets[h$merge[parent,2]]>=dendroZoom$w[1] &&
                    branchCenterOffsets[h$merge[parent,2]]<=dendroZoom$w[2]) {
                    if (dbg>1) cat(' selected\n')
                    selectedClustersIdx<-selectedClustersIdx+1
                    selectedClusters[selectedClustersIdx]<-h$merge[parent,2]
                } else if (dbg>1) cat(' not selected\n')
            }
        }
        if (dbg>1) printVar(selectedClusters)

        selectedClusterCount<-0
        unselectedIndices<-1:df$clusterCount
        for (i in mySeq(1,selectedClustersIdx)) {
            if (dbg>1) printVar(i)
            if (dbg>1) printVar(selectedClusters[i])
            if (!is.na(selectedClusters[i]) && selectedClusters[i]>0) {
                subclusters<-computeSubclusterIndices(h,selectedClusters[i])
                unselectedIndices<-setdiff(unselectedIndices,subclusters)
                if (dbg>1) printVar(subclusters)
                tmp<-clusterId2SegmentIds(subclusters)
                if (dbg>1) printVar(tmp)
                selectedClusterCount<-selectedClusterCount+1
                newClusters[[selectedClusterCount]]$indices<-subclusters
                newClusters[[selectedClusterCount]]$branches<-
                    with(df$allBranches$branches,data.frame(x1s=x1s[tmp],x2s=x2s[tmp],y1s=y1s[tmp],y2s=y2s[tmp]))
            }
        }
        df$clusters<-newClusters[1:selectedClusterCount]
        if (dbg>1) printVar(df$clusters)
        df$unselectedBranches$indices<-unselectedIndices
        idx<-clusterId2SegmentIds(unselectedIndices)
        df$unselectedBranches$branches<-
            with(df$allBranches$branches,data.frame(x1s=x1s[idx],x2s=x2s[idx],y1s=y1s[idx],y2s=y2s[idx]))
    } else {
        selectedClusterCount<-1
        df$clusters<-vector('list',1)
        subclusters<-computeSubclusterIndices(h,df$clusterCount)
        tmp<-clusterId2SegmentIds(subclusters)
        df$clusters[[1]]<-list(indices=subclusters,
            branches=with(df$allBranches$branches,data.frame(x1s=x1s[tmp],x2s=x2s[tmp],y1s=y1s[tmp],y2s=y2s[tmp])))
        tmp<-df$unselectedBranches$indices<-c()
        df$unselectedBranches$branches<-
            with(df$allBranches$branches,data.frame(x1s=x1s[tmp],x2s=x2s[tmp],y1s=y1s[tmp],y2s=y2s[tmp]))
    }

    # compute leaf colors
    df$leafColorIdxs<-computeLeafColorIdxs(df)

    return(list(df=df,selectedClusterCount=selectedClusterCount))
    ### A list of shared data frame 'df' with new cluster selection,
    ### and the number of clusters selected.
}
