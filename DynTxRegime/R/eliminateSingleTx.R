eliminateSingleTx <- function(subsets, ptsSubset){

    useGrps <- which(sapply(X = subsets, FUN = length) > 1L)

    if( length(useGrps) < 0.5 ) {
      UserError("input",
                "No treatment subsets of length > 1.")
    }

    use4fit <- ptsSubset %in% names(subsets)[useGrps]

    if( sum(use4fit) < 0.5 ) {
      UserError("input",
                "No observations with more than 1 tx option.")
    }

    return(use4fit)

}
