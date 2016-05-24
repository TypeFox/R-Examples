print.customizedGlmnet <-
function(x, ...)
{
    cat("\nCall: \n")
    print(x$call)

    cat("\nn =", nrow(x$x$train), "training observations\n")
    cat("m =", nrow(x$x$test), "test observations\n")
    cat("p =", ncol(x$x$train), "predictor variables\n\n")

    for (group in sort(unique(x$groupid))) {
        cat("Model ", group, ": ", sep = "")
        cat(length(x$CTset[[group]]), "train obs. and ")
        cat(sum(x$groupid == group), "test obs.\n")
    }
    cat("\n")
}
