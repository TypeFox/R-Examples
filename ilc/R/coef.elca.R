coef.elca <-
function(object, ...){
    # get the list of parameters
    lcoef <- with(object, list(ax=ax, bx=bx, kt=kt, ag=ag))
    return(structure(lcoef, class='coef'))
}
