# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

### transform canonical vectors back to original scale
#backtransform <- function(A, scale) {
#    apply(A, 2, 
#        function(a) {
#            sa <- a / scale       # divide by scale of corresponding variable
#            sa / sqrt(sum(sa^2))  # divide by norm
#        })
#}

## check if indices are within the limits
checkIndices <- function(indices, max) {
    indices <- as.integer(indices)
    indices[which(indices > 0 & indices <= max)]
}    

## call C++ to compute ranks of observations in a vector (for testing)
fastRank <- function(x) {
    x <- as.numeric(x)
    if(length(x) == 0) return(numeric())        # zero length vector
    .Call("R_rank", R_x=x, PACKAGE="ccaPP")   # call C++ function
}

## get list of control arguments for correlation function
getCorControl <- function(method, control, forceConsistency = TRUE) {
    if(method %in% c("spearman", "kendall", "quadrant")) {
        if(forceConsistency) out <- list(consistent=TRUE)
        else {
            # get default values (the three functions have the same arguments)
            out <- formals(corSpearman)[-(1:2)]
            # check supplied values
            if(is.list(control)) {
                if(!is.null(consistent <- control$consistent)) {
                    out$consistent <- isTRUE(consistent)
                }
            }
        }
    } else if(method == "M") {
        # get default values
        out <- formals(corM)[-(1:2)]
        choices <- eval(out$initial)
        out$initial <- choices[1]
        # check supplied values
        if(is.list(control)) {
            if(!is.null(prob <- control$prob)) {
                out$prob <- as.numeric(prob)
            }
            if(!is.null(initial <- control$initial)) {
                out$initial <- match.arg(initial, choices)
            }
            if(!is.null(tol <- control$tol)) {
                out$tol <- as.numeric(tol)
            }
        }
    } else out <- list()  # this case includes Pearson correlation
    # return list of control arguments
    out
}

## L1 median (for testing)
l1Median <- function(x) {
    # initializations
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if(p == 0) return(numeric())       # no columns
    if(n == 0) return(rep.int(NA, p))  # no observations
    # call C++ function
    .Call("R_l1Median", R_x=x, PACKAGE="ccaPP")
}

### (robustly) standardize the data
#standardize <- function(x, robust = TRUE) {
#    if(robust) {
#        # with median and MAD
#        tmp <- apply(x, 2, function(v) unlist(fastMAD(v)))
#        center <- tmp[1,]  # column medians
#        x <- sweep(x, 2, center, check.margin=FALSE)  # sweep out column centers
#        scale <- tmp[2,]  # column MADs
#        x <- sweep(x, 2, scale, "/", check.margin=FALSE)  # sweep out column scales
#    } else {
#        # with mean and standard deviation
#        center <- colMeans(x)  # compute column means (faster than apply)
#        x <- sweep(x, 2, center, check.margin=FALSE)  # sweep out column centers
#        f <- function(v) sqrt(sum(v^2) / max(1, length(v)-1))
#        scale <- apply(x, 2, f)  # compute column scales with zero means
#        x <- sweep(x, 2, scale, "/", check.margin=FALSE)  # sweep out column scales
#    }
#    # add attributes and return standardized data
#    attr(x, "center") <- center
#    attr(x, "scale") <- scale
#    x
#}
