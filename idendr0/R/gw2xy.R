gw2xy<-function
### Logical grow/width to physical x/y coord conversion.
##keyword<<internal
(
    gw ##<< a list of 'g' and 'w' components
) {
    xy<-gw
    names(xy)<-c('x','y')
    return(xy)
    ### a list of 'x' and 'y' components
}
