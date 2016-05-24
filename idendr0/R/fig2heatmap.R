fig2heatmap<-function
### Figure -> heatmap coord conversion.
##keyword<<internal
(
    gw ##<< a list of 'g' and 'w' components
) {
    dbg.tx<-.gfc(dbg.tx,required=FALSE)
    if (!is.null(dbg.tx) && dbg.tx) print('fig2heatmap called')

    return(gw)
    ### a list of 'g' and 'w' components
}
