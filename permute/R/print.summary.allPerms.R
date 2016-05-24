`print.summary.allPerms` <- function(x, ...) {
    dims <- dim(x)
    control <- attr(x, "control")
    observed <- attr(x, "observed")
    attributes(x) <- NULL
    dim(x) <- dims
    cat("\n")
    writeLines(strwrap("Complete enumeration of permutations\n",
        prefix = "\t"))
    print(control)
    cat("\nAll permutations:\n")
    writeLines(paste("Contains observed ordering?:", ifelse(observed, "Yes", "No"),
              "\n"))
    print(x)
    return(invisible(x))
}
