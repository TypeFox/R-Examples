#'@S3method print bootnodes
print.bootnodes <-
function(x, ...)
{
    cat("BOOTSTRAP FOR PATH COEFFICIENTS OF TERMINAL NODES", "\n")
    cat("\n")
    cat("--------------------------------------------------------", "\n") 
    cat("$PC: ORIGINAL PATH COEFFICIENTS", "\n\n")
    print(x$PC, digits=3, print.gap=2)
    cat("\n")
    cat("--------------------------------------------------------", "\n") 
    cat("$PMB: MEAN VALUE FOR BOOTSTRAP PATH COEFFICIENTS", "\n\n")
    print(x$PMB, digits=3, print.gap=2)
    cat("\n")
    cat("--------------------------------------------------------", "\n") 
    cat("$PSB: STANDARD ERROR FOR BOOTSTRAP PATH COEFFICIENTS", "\n\n")
    print(x$PSB, digits=3, print.gap=2)
    cat("\n")
    cat("--------------------------------------------------------", "\n") 
    cat("$PP05: PERCENTILE 5 FOR BOOTSTRAP PATH COEFFICIENTS", "\n\n")
    print(x$PP05, digits=3, print.gap=2)
    cat("\n")
    cat("--------------------------------------------------------", "\n") 
    cat("$PP95: PERCENTILE 95 FOR BOOTSTRAP PATH COEFFICIENTS", "\n\n")
    print(x$PP95, digits=3, print.gap=2)
    invisible(x)
}

