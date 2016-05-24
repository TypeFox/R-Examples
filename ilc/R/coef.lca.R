coef.lca <-
function(object,...){
    # get the list of parameters
    lcoef <- with(object, list(ax=ax, bx=bx, kt=kt))
    return(structure(lcoef, class='coef'))
}
