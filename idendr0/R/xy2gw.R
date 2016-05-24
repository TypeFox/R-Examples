xy2gw<-function
### Physical x/y to logical grow/width coord conversion.
##keyword<<internal
(
    xy ##<< a list of 'x' and 'y' components

) {
    gw<-xy
    names(gw)<-c('g','w')
    gw
    ### a list of 'g' and 'w' components
}
