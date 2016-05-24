print.metabdata <- function(x, ...)
{
    x$output <- editcolnames(x$output)
    if (dim(x$output)[2] > 5) {
        prettydata <- cbind(
            head(x$output)[1:10], rep('...',dim(head(x$output))[1])
        )
        colnames(prettydata)[11] <- '...'
        print (prettydata)
        cat("\nOnly the first few rows and columns are printed.\n")
    } else {
        print(head(x$output))
        cat("\nOnly the first few rows are printed.\n")
    }
    cat("\nGroups:\n")
    print(x$groups)
    cat("\nSamples:\n")
    print(x$samples)
}
