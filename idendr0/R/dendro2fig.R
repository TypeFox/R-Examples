dendro2fig<-function
### Dendrogram -> figure coords conversion.
##keyword<<internal
(
    gw ##<< a list of 'g' and 'w' components, g=gw[1]: dendro height ... 0,
    ## w=gw[2]: n ... 1 (charm to strange)
) {
    dbg.tx<-.gfc(dbg.tx,required=FALSE)
    if (!is.null(dbg.tx) && dbg.tx) print('dendro2fig called')

    names(gw)<-c('g','w')
    return(gw)
    ### converted values, a list of 'g' and 'w' components
}
