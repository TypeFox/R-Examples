`fevaluation` <-
function(z,ListOfFunctions, useApply = FALSE){
    # This is the evaluation:
    f <- function(f1, z){ apply(z, 2, f1) }
    tF <- matrix(unlist(lapply(ListOfFunctions, f, z = z)), nrow = length(z[1,]))
    F <- t(tF)
    F
}

