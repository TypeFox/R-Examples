## Modification of function Anova in package genefilter
oneWayAnova <- function(cov, na.rm = TRUE, var.equal = FALSE){
    function(x) {
        if (na.rm) {
            drop <- is.na(x)
            x <- x[!drop]
            cov <- cov[!drop]
        }
        oneway.test(x ~ cov, var.equal = var.equal)$p.value
    }
}
