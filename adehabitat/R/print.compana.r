"print.compana" <- function(x, ...)
{
    if (!inherits(x, "compana"))
        stop("should be an object of class \"compana\"")
    cat("************ Compositional analysis of habitat use ***************\n\n")
    cat("The analysis was carried out with", nrow(x$used),
        "animals and", ncol(x$used), "habitat types\n")
    cat("1. Test of the habitat selection:\n")
    cat("  ", x$type.test, "test\n")
    print(x$test)
    cat("\n2. Ranking of habitats (profile):\n")
    print(x$profile, quote=FALSE)
}

