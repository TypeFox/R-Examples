clusterId2SegmentIds<-function
### Convert cluster IDs to indices of dendrogram segments.
##keyword<<internal
(
    ids ##<< cluster IDs
) {
    return(rep(3*(ids-1),each=3)+(1:3))
    ### segment IDs
}
