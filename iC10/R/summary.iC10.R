summary.iC10 <- function(object,...) {
    cat("\nDistribution of samples\n")
    print(table(object$class))
    cat("\nSummary of probablities of assignment for all samples:\n")
    tmp <- apply(object$posterior, 2, max)
    summary(tmp)
}
