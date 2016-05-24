computeLeafOrder<-function
### Compute the assignment of dendrogram leafs to the original
### observations.
##keyword<<internal
(
    merging, ##<< (n-1) by 2 matrix describing HCA merging, usually the
    ## 'merge' component of the return value from 'hclust'
    dbg=0  ##<< debug verbosity level
) {
    # m = # of clusters
    m<-dim(merging)[1]
    if (dbg) printVar(m)
    # n = # of leafs
    n<-m+1
    if (dbg) printVar(m)

    # compute 'order'
    lifo<-rep(0,m)
    lifoIdx<-0
    ordering<-rep(0,n)
    # running leaf index being assigned to leafs in order
    oi<-1
    # push the top merging
    if (dbg) cat(sprintf('pushing %d\n',m))
    lifoIdx<-lifoIdx+1
    lifo[lifoIdx]<-m
    while (lifoIdx>0) {
        clstr<-lifo[lifoIdx]
        lifoIdx<-lifoIdx-1
        if (dbg) cat(sprintf('popping %d\n',clstr))
        if (clstr<0) {
            # leaf found -> assign it an index
            if (dbg) cat(sprintf('  -> %d\n',oi))
            ordering[-clstr]<-oi
            oi<-oi+1
        } else {
            # process both subclusters
            for (i in 2:1) {
                subclstr<-merging[clstr,i]
                if (dbg) cat(sprintf('pushing %d\n',subclstr))
                lifoIdx<-lifoIdx+1
                lifo[lifoIdx]<-subclstr
            }
        }
    }
    if (dbg) printVar(ordering)
    return(ordering)
    ### leaf ordering
}
