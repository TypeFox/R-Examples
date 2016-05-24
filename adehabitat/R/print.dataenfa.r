"print.dataenfa" <- function (x, ...)
{
    if (!inherits(x, "dataenfa"))
        stop("Object of class 'dataenfa' expected")
    cat("Data ENFA\n")
    cat("\n List of 4 elements:\n\n")
    sumry <- array("", c(1, 4), list(1, c("data.frame", "nrow", "ncol",
        "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "table of pixels")
    class(sumry) <- "table"
    print(sumry)
    cat("\n")
    sumry <- array("", c(2, 3), list(1:2, c("vector", "length", "content")))
    sumry[1, ] <- c("$pr", length(x$pr), "vector of presence")
    sumry[2, ] <- c("$index", length(x$index), "position of the rows")
    class(sumry) <- "table"
    print(sumry)
    cat("\n$attr: attributes of the initial kasc\n")
}

