smoothHeatmap<-function
### Smooth heatmap by averaging data over observations in clusters.
##keyword<<internal
(
    x, ##<< a data frame holding observations tha were clustered

    ch, ##<< number of clusters specifying the amount of smoothing:
    ## observations in clusters will get smoothed together, the value
    ## of 'n' specifies no smoothing while the value of 1 would lead
    ## to maximal smoothing

    dbg.heatmap = 0 ##<< debug verbosity level
) {
    
    i0<-0
    if (dbg.heatmap) chProcessed<-rep(FALSE,max(ch))
    for (i in 1:max(ch)) {
        if (dbg.heatmap>2) printVar(i)
        if (dbg.heatmap>2) printVar(ch[i0+1])
        if (dbg.heatmap) if (chProcessed[ch[i0+1]]) {cat('error: invalid HCA \'order\'');}#browser()}
        if (dbg.heatmap) chProcessed[ch[i0+1]]<-TRUE
        cnt<-sum(ch==ch[i0+1])
        if (dbg.heatmap>2) printVar(i0)
        if (dbg.heatmap>2) printVar(cnt)
        if (dbg.heatmap>2) printVar(i0+(1:cnt))
        if (dbg.heatmap>2) printVar(x[i0+(1:cnt),])
        x[i0+(1:cnt),]<-rep(colMeans(x[i0+(1:cnt),,drop=FALSE]),each=cnt)
        i0<-i0+cnt
    }
    if (dbg.heatmap>2) printVar(x)
    return(x)
    ### a data frame holding smoothed observations
}
