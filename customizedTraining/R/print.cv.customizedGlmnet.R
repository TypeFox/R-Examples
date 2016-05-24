print.cv.customizedGlmnet <-
function(x, ...)
{
    cat("\nCall: \n")
    print(x$call)

    cat("\nn =", nrow(x$fit$x$train), "training observations\n")
    cat("m =", nrow(x$fit$x$test), "test observations\n")
    cat("p =", ncol(x$fit$x$train), "predictor variables\n")

    cat("\nNumber of groups chosen:", x$G.min, "\n\n")

    for (group in sort(unique(x$fit$groupid))) {
        cat("Model ", group, ": ", sep = "")
        cat(length(x$fit$CTset[[group]]), "train obs. and ")
        cat(sum(x$fit$groupid == group), "test obs.\n")

        cat(" No. of variables selected: ")
        if (is.element(x$fit$family, c("gaussian", "binomial"))) {
            cat(length(x$selected[[group]][[1]]))
        } else if (x$fit$family == "multinomial") {
            for (level in levels(x$fit$y)) {
                cat("Class ", level, ": ", sep = "")
                cat(length(x$selected[[group]][[level]][[1]]))
                cat("   ")
            }
        }
        cat("\n\n")
    }
}
