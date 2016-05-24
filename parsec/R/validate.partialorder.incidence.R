validate.partialorder.incidence <-
function(m)  {
    warn <- NULL
    #if(!is.partialorder(m)) stop("not an incidence martix")
    if(!binary(m)) warn <- c(warn, "binary\n")
    if(!antisymmetry(m)) warn <- c(warn, "antisymmetric\n")
    if(!reflexivity(m)) warn <- c(warn, "refelxive\n")
    if(!transitivity(m)) warn <- c(warn, "transitive\n")
    if(!is.null(warn)) stop("the matrix is not:\n", warn)
    class(m) <- "incidence"
    return(m)
}
