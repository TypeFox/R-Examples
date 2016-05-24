heatmap2fig<-function
### Heatmap -> figure coord conversion.
##keyword<<internal
(
    gw ##<< a list of 'g' and 'w' components
) {
    dbg.tx<-.gfc(dbg.tx,required=FALSE)
    if (!is.null(dbg.tx) && dbg.tx) print('heatmap2fig called')

    return(gw)
    ### a list of 'g' and 'w' components
}
