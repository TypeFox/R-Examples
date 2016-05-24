coef.rh <-
function(object, ...){
    # get the list of parameters
    lcoef <- with(object, list(ax=ax, bx0=bx0, bx1=bx, kt=kt, itx=itx))
    # check for parameters that were excluded from the fitting
    ind <- NULL
    for(i in seq(lcoef)) if(!any(bool(lcoef[[i]]))) ind <- c(ind, i)
    if (!is.null(ind)) lcoef <- lcoef[-ind]
    return(structure(lcoef, class='coef'))
}
