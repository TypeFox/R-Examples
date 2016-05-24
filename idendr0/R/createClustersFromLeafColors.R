createClustersFromLeafColors<-function
### Initialize clusters (i.e. cluster selection) from leaf colors (i.e.
### assignment of some observations to some clusters)
##keyword<<internal
(
    df, ##<< shared data frame, see 'prepareDendro'

    leafColors, ##<< leaf color indices, i.e. a nonnegative number
    ## assigned to each observation. 0 means given observation is not
    ## member of any cluster, a value i>0 means given observation is a
    ## member of cluster 'i'.

    maxClusterCount, ##<< max cluster count

    dbg=0 ##<< debug verbosity level
    ) {

    if (dbg) cat('createClustersFromLeafColors called\n')
    if (dbg>1) printVar(leafColors)

    # make commonly used variables accessible more easily
    n<-df$n
    clusterCount<-df$clusterCount
    h<-df$h

    if (!is.null(leafColors) && length(leafColors)==n &&
        is.numeric(leafColors) && any(leafColors>0)) {
        # initialize clusters from user-supplied initial assignment
        if (dbg) cat('Initializing clusters from user-supplied clusters argument...\n')
        if (any(leafColors<0)) {
            stop("Invalid initial clusters: numeric vector holding nonnegative numbers needed.")
        }
        df$leafColorIdxs<-leafColors
        # assignment of each leaf and subcluster to user-supplied cluster IDs
        clst<-c(leafColors,rep(0,clusterCount))
        # propagate clusters assignments up the hierarchy
        mergeMember2Idx<-function(m) ifelse(m<0,-m,n+m)
        for (i in 1:clusterCount) {
            if (clst[mergeMember2Idx(h$merge[i,1])]==clst[mergeMember2Idx(h$merge[i,2])]) {
                # color cluster according to the color of both its subclusters
                clst[n+i]<-clst[mergeMember2Idx(h$merge[i,1])]
            }
        }
        if (dbg>1) printVar(clst)

        # traverse the hierarchy from top and check that user-supplied
        # colors match real clusters

        # indices of already processed cluster IDs
        clusterIdsProcessed<-c()
        # indices of supplied cluster IDs
        clustersIdsSupplied<-sort(setdiff(leafColors,0))
        if (dbg) printVar(clustersIdsSupplied)
        if (max(clustersIdsSupplied)>maxClusterCount) {
            stop(sprintf('user-supplied cluster %d greater than maxClusterCount %d.',
                max(clustersIdsSupplied),maxClusterCount))
        }
        unselectedIndices<-1:n
        clst<-clst[n+(1:clusterCount)]
        for (i in seq(clusterCount,1,-1)) {
            # is this cluster assigned to some user-supplied cluster
            # (and has not been processed yet)?
            if (clst[i]!=0 && !(clst[i] %in% clusterIdsProcessed)) {
                if (dbg) printVar(clst[i])
                # does this cluster contain exactly the leafs as
                # specified by user?
                leafs<-computeMemberIndices(h,i)
                if (dbg>1) printVar(leafs)
                if (dbg>1) printVar(which(leafColors[leafs]==clst[i]))
                if (length(leafs)==sum(leafColors==clst[i]) && all(leafColors[leafs]==clst[i])) {
                    if (dbg) cat('cluster',clst[i],'is OK\n')
                    subclusters<-computeSubclusterIndices(h,i)
                    unselectedIndices<-setdiff(unselectedIndices,subclusters)
                    tmp<-clusterId2SegmentIds(subclusters)
                    df$clusters[[clst[i]]]<-list(indices=subclusters,
                        branches=with(df$allBranches$branches,data.frame(x1s=x1s[tmp],x2s=x2s[tmp],y1s=y1s[tmp],y2s=y2s[tmp])))
                } else {
                    stop(paste("Invalid 'clusters' argument: observations of ID",clst[i],"do not form a cluster."))
                }
                clusterIdsProcessed<-c(clusterIdsProcessed,clst[i])
                if (length(clusterIdsProcessed)==length(clustersIdsSupplied)) {
                    # we're done: no more clusters to process
                    break
                }
            }
        }
        # mark remaining (unselected) indices/branches as unselected
        tmp<-clusterId2SegmentIds(unselectedIndices)
        df$unselectedBranches<-list(indices=unselectedIndices,
            branches=with(df$allBranches$branches,data.frame(x1s=x1s[tmp],x2s=x2s[tmp],y1s=y1s[tmp],y2s=y2s[tmp])))
        df$leafColorIdxs<-leafColors
    } else {
        if (dbg) cat('User-supplied leaf colors do not contain any clusters.\n')
    }

    return(df)
    ### shared data frame 'df' with cluster selection based on 'leafColors'
}
