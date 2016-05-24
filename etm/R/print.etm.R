print.etm <- function(x, covariance = TRUE, whole = TRUE, ...) {
    if (!inherits(x, "etm"))
        stop("'x' must be of class 'etm'")
    absorb <- setdiff(levels(x$trans$to), levels(x$trans$from))
    transient <- unique(x$state.names[!(x$state.names %in% absorb)])
    cat(paste("Multistate model with", length(transient), "transient state(s)\n",
              "and", length(absorb), "absorbing state(s)\n\n", sep = " "))
    cat("Possible transitions:\n")
    print(x$trans, row.names = FALSE)
    cat("\n")
    cat(paste("Estimate of P(", x$s, ", ", x$t, ")\n", sep = ""))
    print(x$est[, , dim(x$est)[3]]); cat("\n")
    if (!is.null(x$cov) & covariance == TRUE) {
        if (whole) {
            cat(paste("Estimate of cov(P(", x$s, ", ", x$t, "))\n", sep = ""))
            print(x$cov[, , dim(x$cov)[3]])
        }
        else {
            cov <- x$cov[, , dim(x$cov)[3]][rowSums(x$cov[, , dim(x$cov)[3]]) != 0, ]
            cova <- cov[, colSums(cov) != 0]
            cat(paste("Estimate of cov(P(", x$s, ", ", x$t, "))\n", sep = ""))
            print(cova)
        }
    }
    invisible()
}
    
    
            
