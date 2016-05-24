## as.matrix.permutationMatrix - an S3 method to convert to the S3
## matrix class. Essentially this just strips attributes and updates
## the class to only "matrix"

`as.matrix.permutationMatrix` <- function(x, ...) {
    attr(x, "seed") <- NULL
    attr(x, "control") <- NULL
    class(x) <- "matrix"
    x
}
