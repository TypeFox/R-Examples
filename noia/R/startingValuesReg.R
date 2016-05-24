startingValuesReg <-
function (reg, max.level = 2, max.dom = 2, e.unique = FALSE, 
    nloc = NULL) 
{
    if (is.null(reg)) {
        return(startingValuesNothing(nloc, max.level, max.dom, 
            e.unique))
    }
    else if (class(reg) == "noia.linear") {
        return(startingValuesLinear(reg, max.level, max.dom, 
            e.unique))
    }
    else if (class(reg) == "noia.multilinear") {
        return(startingValuesMultilinear(reg, max.level, max.dom, 
            e.unique))
    }
    else {
        stop("Object of class \"noia.linear\" or \"noia.multilinear\" expected\n")
    }
}
