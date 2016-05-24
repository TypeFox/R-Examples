ghq <- function(n.points = 1, modified = TRUE){
    weights <- numeric(n.points)
    zeros <- numeric(n.points)
    res <- .Fortran("ghq",
                    as.integer(n.points),
                    zeros = as.double(zeros),
                    weights = as.double(weights),
                    as.logical(modified),
                    PACKAGE = "eha")
    list(weights = res$weights, zeros = res$zeros)
}
