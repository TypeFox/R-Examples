edges <-
function(plants, width) {
# Expand or contract window, to handle edge effects
    if(width < 0)
        return(plants[erosion(as.owin(plants), -width)])
    if(width > 0)
        return(periodify(plants)[dilation(as.owin(plants), width)])
    return(plants)
}
