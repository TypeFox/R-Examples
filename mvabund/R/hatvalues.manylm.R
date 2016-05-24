hatvalues.manylm <- 
function (model, infl = manylm.influence(model, do.coef = FALSE), 
    ...) {

    hat <- infl$hat
    if(is.null(names(hat)))
    	names(hat) <- rownames(infl$wt.res)
    return(hat)

}
