c.WGassociation<-
function (...)
{
    allargs <- list(...)
    allargs <- allargs[sapply(allargs, length) > 0]
    n <- length(allargs)
    if (n == 0)
        return(structure(list(), class = "data.frame", row.names = integer()))
    lapply(1:n, function(i) if (!inherits(allargs[[i]], "WGassociation"))
        stop("Please supply 'WGassociation' objects"))
    x1 <- allargs[[1]]
    for (i in 2:n) {
        xi <- allargs[[i]]
        if (any(attr(x1, "models") != attr(xi, "models")))
            stop("All objects should have identical structure")
        if (attr(x1, "quantitative") != attr(xi, "quantitative"))
            stop("All objects should have identical structure")
    }
    for (i in 2:n) {
        xi <- allargs[[i]]
        out <- rbind(attr(x1, "pvalues"), attr(xi, "pvalues"))
        attr(out, "tables") <- c(attr(x1, "tables"), attr(xi,
            "tables"))
        attr(out, "label.SNPs") <- c(attr(x1, "label.SNPs"),
            attr(xi, "label.SNPs"))
        attr(out, "colSNPs") <- c(attr(x1, "colSNPs"), attr(xi,
            "colSNPs"))
        attr(out, "gen.info") <- rbind(attr(x1, "gen.info"),
            attr(xi, "gen.info"))
        attr(out, "pvalues") <- out
        x1<-out
    }
    class(out) <- c("WGassociation", "data.frame")
    out
}