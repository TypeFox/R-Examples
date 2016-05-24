prepareDendro<-function
### Perform all computations needed to display a dendrogram and
### initialize a shared data frame.
###
### This is an internal function not to be called by the user.
##keyword<<internal
(
    h, ##<< an object of class 'stats::hclust' describing a clustering

    x = NULL, ##<< a data frame holding observations tha were clustered
    ## giving rise to 'h', converted to numeric, and scaled

    xOrig = NULL, ##<< a data frame holding observations tha were clustered
    ## giving rise to 'h'

    doFlipG = TRUE, ##<< should branches' heights ("grow") be flipped
    ## such that the heights of elementary observations correspond not
    ## to 0, but to the height of the dendrogram?

    dbg = 0 ##<< debug verbosity level
) {

    df<-NULL
    df$doFlipG<-doFlipG
    df$h<-h
    df$x<-x
    df$clusterCount<-clusterCount<-length(h$height)
    df$n<-n<-df$clusterCount+1
    if (is.null(x)) {
        df$k<-0
    } else {
        df$k<-ncol(x)
    }

    if (dbg) printVar(n)
    if (dbg) printVar(clusterCount)

    if (dbg) cat('Computing memberCount...\n')
    # memberCount holds the number of elements in clusters, it
    # spans the range of elements in h$merge ranging from -n to n-1,
    # thus:
    # - memberCount [1..n] (corresponding to -n..-1)
    #       - holds 1 (trivial clusters of size 1)
    # - memberCount [n+1] (corresponding to 0)
    #       - holds NA (undefined, there is no cluster of ID 0)
    # - memberCount [n+1+(1..clusterCount)] (corresponding to
    #       1..clusterCount 0)
    #       - holds the number of observations in non-trivial clusters
    memberCount<-c(rep(1,n),NA,rep(NA,clusterCount))
    for (i in 1:clusterCount) {
        memberCount[n+1+i]<-memberCount[n+1+h$merge[i,1]]+memberCount[n+1+h$merge[i,2]]
    }
    if (dbg>1) printVar(memberCount)

    if (dbg) cat('Computing prototypes...\n')
    # prototypes are single elementary observations selected from each
    # cluster; the purpose of prototypes is to ease determinating the
    # order of clusters in `h$order'
    prototypes<-c(n:1,NA,rep(NA,clusterCount))
    for (i in 1:clusterCount) {
        prototypes[n+1+i]<-prototypes[n+1+h$merge[i,1]]
    }
    if (dbg>1) printVar(prototypes)

    # compute the order of observations in h$order
    # (used to determine which of two observations comes first in
    # h$order)
    orderOfOrder<-rep(NA,n)
    orderOfOrder[h$order]<-1:n

    if (dbg) cat('Computing offsets...\n')
    clusterOffsets<-rep(NA,clusterCount)
    clusterOffsets[clusterCount]<-0
    leafOffsets<-rep(NA,n)
    for (i in clusterCount:1) {
        # determine whether members of h$merge[i,1] or h$merge[i,2]
        # come first in h$order
        if (orderOfOrder[prototypes[n+1+h$merge[i,1]]] < orderOfOrder[prototypes[n+1+h$merge[i,2]]]) {
            # stack h$merge[i,1] first
            c1<-h$merge[i,1]
            c2<-h$merge[i,2]
        } else {
            # stack h$merge[i,2] first
            c1<-h$merge[i,2]
            c2<-h$merge[i,1]
        }
        # lower branch offset is dictated by the parent branch
        if (c1>0) {
            clusterOffsets[c1]<-clusterOffsets[i]
        } else {
            leafOffsets[n+1+c1]<-clusterOffsets[i]
        }
        # upper branch offset is dictated by the sum of parent branch
        # offset and the width of the lower branch
        if (c2>0) {
            clusterOffsets[c2]<-clusterOffsets[i]+memberCount[n+1+c1]
        } else {
            leafOffsets[n+1+c2]<-clusterOffsets[i]+memberCount[n+1+c1]
        }
    }
    rm(memberCount)
    if (dbg>1) printVar(clusterOffsets)
    rm(clusterOffsets)
    if (dbg>1) printVar(leafOffsets)

    if (dbg) cat('Computing yPos...\n')
    yPos<-c(1+leafOffsets,NA,rep(NA,clusterCount))
    rm(leafOffsets)
    for (i in 1:clusterCount) {
        yPos[n+1+i]<-mean(yPos[n+1+h$merge[i,]])
    }
    if (dbg>1) printVar(yPos)

    #leafOrder<-order(computeLeafOrder(h$merge))
    leafOrder<-h$order

    if (dbg) cat('Computing segments...\n')
    x1s<-x2s<-y1s<-y2s<-rep(NA,3*clusterCount)
    ii<-0
    branchCenterOffsets<-rep(NA,clusterCount)
    # heights of both trivial clusters of size 1 and non-trivial clusters
    heights<-c(rep(0,n),NA,h$height)
    topClusterHeight<-h$height[clusterCount]
    for (i in 1:clusterCount) {
        x0<-h$height[i]
        x1<-heights[n+1+h$merge[i,1]]
        x2<-heights[n+1+h$merge[i,2]]
        if (doFlipG) {
          # mirror Xs: x=0 corresponding to the top-most cluster, x>0 to leafs
          x0<-topClusterHeight-x0
          x1<-topClusterHeight-x1
          x2<-topClusterHeight-x2
        }

        y1<-yPos[n+1+h$merge[i,1]]
        y2<-yPos[n+1+h$merge[i,2]]
        branchCenterOffsets[i]<-mean(c(y1,y2))
        # x1,y1 -> x0,y1
        # x0,y1 -> x0,y2
        # x2,y2 -> x0,y2
        x1s[ii+(1:3)]<-c(x1,x0,x2)
        x2s[ii+(1:3)]<-c(x0,x0,x0)
        y1s[ii+(1:3)]<-c(y1,y1,y2)
        y2s[ii+(1:3)]<-c(y1,y2,y2)
        ii<-ii+3
    }
    rm(heights)
    if (dbg>1) printVar(x1s)
    if (dbg>1) printVar(x2s)
    if (dbg>1) printVar(y1s)
    if (dbg>1) printVar(y2s)
    #coords1<-gw2xy(dendro2fig(with(df$unselectedBranches$branches,list(x1s,y1s))))
    #coords2<-gw2xy(dendro2fig(with(df$unselectedBranches$branches,list(x2s,y2s))))
    coords1<-gw2xy(dendro2fig(list(x1s,y1s)))
    coords2<-gw2xy(dendro2fig(list(x2s,y2s)))
    x1s<-coords1[[1]]
    y1s<-coords1[[2]]
    x2s<-coords2[[1]]
    y2s<-coords2[[2]]
    if (dbg>1) printVar(x1s)
    if (dbg>1) printVar(x2s)
    if (dbg>1) printVar(y1s)
    if (dbg>1) printVar(y2s)
    branchCenterOffsets<-dendro2fig(list(NULL,branchCenterOffsets))[[2]]

    if (dbg) cat('Initializing the shared data frame...\n')
    df$unselectedBranches<-df$allBranches<-list(indices=1:clusterCount,branches=data.frame(x1s=x1s,x2s=x2s,y1s=y1s,y2s=y2s))
    df$clusters<-NULL
    df$currentCluster<-1
    df$clusterCount<-clusterCount
    df$branchCenterOffsets<-branchCenterOffsets
    df$leafOrder<-leafOrder
    df$leafColorIdxs<-rep(0,n)
    df$xOrdered<-x[leafOrder,,drop=FALSE]
    df$xOrderedSmoothed<-df$xOrdered
    if (is.null(xOrig)) {
        df$xOrigOrdered<-NULL
    } else {
        df$xOrigOrdered<-xOrig[leafOrder,,drop=FALSE]
    }
    df$elemClusterCount<-df$n

    return(df)
    ### shared data frame
}
