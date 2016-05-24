fig2dendro<-function
### Fig -> dendro coord conversion.
##keyword<<internal
(
    gw ##<< ##<< a list of 'g' and 'w' components
) {
    dbg.tx<-.gfc(dbg.tx,required=FALSE)
    if (!is.null(dbg.tx) && dbg.tx) print('fig2dendro called')

    return(gw)
    ### converted values, a list of 'g' and 'w' components
}
