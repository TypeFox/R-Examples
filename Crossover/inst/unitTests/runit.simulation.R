simulation <- function (result, refDesign, model) {
# There's probably some examples, but there are some examples of people using
#   solve(t(X) %*% W %*% X) %*% W %*% Y
# to compute regression coefficients, too.
#   -- Thomas Lumley (discussing usefulness of evaluation order in lapply) - R-help (March 2006)
}

test.simulation <- function() {
    if (!"extended" %in% strsplit(Sys.getenv("CROSSOVER_UNIT_TESTS"),",")[[1]]) {
        cat("Skipping simulation.\n")
        return()
    }
    # Example
    mu <- c(0.1, 0.2, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06)
    v <- 4
    
    model <- 1
    
}

