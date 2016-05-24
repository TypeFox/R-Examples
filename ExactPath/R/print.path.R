print.path = function(x, ...)
{
    cat("breaks:\n")
    print(round(x$breaks, 4))
    cat("\nIndicator matrix (parameters x breaks):\n")
    print(x$tau)
    cat("\nBeta matrix (parameters x breaks):\n")
    print(round(x$beta, 4))
    cat("\nScore matrix (parameters x breaks):\n")
    print(round(x$score, 4))
}
