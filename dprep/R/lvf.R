lvf <-
function (data, lambda, maxiter) 
{
    f <- dim(data)[2] - 1
    bestsize <- f
    variables <- 1:f
    bestsubset <- 1:f
    count <- 0
    bestinselec = lambda
    for (i in 1:maxiter) {
        repeat {
            indic <- rbinom(f, 1, 0.5)
            if (sum(indic) > 0) 
                break
        }
        subset <- variables[indic > 0]
        csubset <- length(subset)
        if (csubset <= bestsize) {
            inselec <- inconsist(data[, c(subset, f + 1)])
            if (inselec < bestinselec) {
                if (csubset < bestsize) {
                  bestsubset <- subset
                  bestsize <- csubset
                }
            }
        }
    }
    cat("The inconsistency of the best subset is\n")
    cat(inconsist(data[, c(bestsubset, f + 1)]))
    cat("\nThe best subset of features is:\n")
    bestsubset
}
