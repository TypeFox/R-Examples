# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

inclusionProb <- function(prob, size) {
    prob <- as.numeric(prob)
    size <- as.integer(size[1])
    if(length(prob) == 0) return(numeric())
    if(length(size) == 0 || size == 0) return(rep.int(0, length(prob)))
    if(size < 0) stop("'size' must be a non-negative integer")
#    .Call("inclusionProb", prob, size, PACKAGE="simFrame")
    .Call("R_inclusionProb", R_prob=prob, R_size=size, PACKAGE="simFrame")
}
