print.noia.vardec <-
function (x, ...) 
{
    cat("\nVariance decomposition:\n")
    if (!is.null(x$V_G)) {
        totvar <- x$V_G
    }
    else {
        totvar <- sum(unlist(x), na.rm = TRUE)
    }
    cat(paste("\tTotal genetic variance:", format(totvar, digits = 5), 
        "\n"))
    for (level in setdiff(names(x), "V_G")) {
        cat(paste("\tOrder", level, "\tTotal:\t", format(sum(x[[level]], 
            na.rm = TRUE), digits = 5, nsmall = 4), "\t(", format(100 * 
            sum(x[[level]], na.rm = TRUE)/totvar, digits = 3, 
            nsmall = 1), "%)", "\n"))
        for (i in 1:(length(x[[level]]))) {
            if (is.finite(x[[level]][i])) {
                cat(paste("\t\t\t", names(x[[level]])[i], "\t", 
                  format(x[[level]][i], digits = 5, nsmall = 4), 
                  "\t(", format(100 * x[[level]][i]/totvar, digits = 3, 
                    nsmall = 1), "%)", "\n"))
            }
        }
    }
}
