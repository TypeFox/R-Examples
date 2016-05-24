validation.nor <-
function (no.nor, mean.vec.nor = NULL, var.nor = NULL) 
{
    if ((no.nor < 0) | (floor(no.nor) != no.nor)) {
        stop("Number of normal variables \nmust be an integer whose value or 0 !\n")
    }
    if (!is.null(mean.vec.nor) & no.nor == 0) {
        stop("Mean vector for the normal part is specified while no.nor=0!\n")
    }
    if (!is.null(var.nor) & no.nor == 0) {
        stop("Vector of variances for the normal part is specified while no.nor=0!\n")
    }
    if (is.null(mean.vec.nor) & no.nor > 0) {
        stop("Mean vector for the normal part is not specified while no.nor>0!\n")
    }
    if (is.null(var.nor) & no.nor > 0) {
        stop("Vector of variances for the normal part is not specified while no.nor>0!\n")
    }
    if (!is.null(mean.vec.nor) & !is.null(var.nor) & no.nor > 
        0) {
        if (length(mean.vec.nor) != no.nor) {
            stop("Mean vector for the normal part is misspecified, \ndimension is wrong!\n")
        }
        if (length(var.nor) != no.nor) {
            stop("Vector of variances for the normal part is misspecified, \ndimension is wrong!\n")
        }
        if (min(var.nor) <= 0) {
            stop("Variances must be positive!\n")
        }
    }
cat("No problems are detected for the marginal specification of normal variables! \n")
}

