validation.bin <-
function (no.bin, prop.vec.bin = NULL) 
{
    if ((no.bin < 0) | (floor(no.bin) != no.bin)) {
        stop("Number of binary variables \nmust be a non-negative integer\n")
    }
    else if (!is.null(prop.vec.bin)) {
        if (no.bin == 0) {
            stop("Proportion vector is specified while no.bin=0")
        }
        else if ((min(prop.vec.bin) <= 0) | (max(prop.vec.bin) >= 
            1)) {
            stop("Proportions for binary variables must be between 0 and 1!\n")
        }
        else if (length(prop.vec.bin) != no.bin) {
            stop("Proportion vector is misspecified, dimension is wrong!\n")
        }
    }
    else if (is.null(prop.vec.bin)) {
        if (no.bin > 0) {
            stop("Proportion vector is not specified while no.bin > 0")
        }
    }
cat("No problems are detected for the marginal specification of binary variables! \n")
}

