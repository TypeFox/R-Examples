print.snp<-
function (x, quote = FALSE, max.levels = NULL, width = getOption("width"),
    ...)
{
    ord <- is.ordered(x)
    if (length(x) <= 0)
        cat(if (ord)
            "ordered"
        else "factor", "(0)\n", sep = "")
    else {
        xx <- x
#        class(xx) <- NULL
#       levels(xx) <- NULL
        attributes(xx)<-NULL
        xx[] <- as.character(x)
        print(xx, quote = quote, ...)
    }
    maxl <- if (is.null(max.levels))
        TRUE
    else max.levels
    if (maxl) {
        n <- length(lev <- encodeString(levels(x), quote = ifelse(quote,
            "\"", "")))
        colsep <- if (ord)
            " < "
        else " "
        T0 <- "Genotypes: "
        if (is.logical(maxl))
            maxl <- {
                width <- width - (nchar(T0, type = "w") + 3 +
                  1 + 3)
                lenl <- cumsum(nchar(lev, type = "w") + nchar(colsep,
                  type = "w"))
                if (n <= 1 || lenl[n] <= width)
                  n
                else max(1, which(lenl > width)[1] - 1)
            }
        drop <- n > maxl
        cat(if (drop)
            paste(format(n), ""), T0, paste(if (drop)
            c(lev[1:max(1, maxl - 1)], "...", if (maxl > 1) lev[n])
        else lev, collapse = colsep), "\n", sep = "")
        cat("Alleles: ",attr(x,"allele.names"),"\n")
    }
    invisible(x)
}
