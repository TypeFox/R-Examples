.apermAssay <- function(TableChar, fun, log,
                        inversPerm, dimension, names,
                        combinedTreatment = FALSE,
                        Inner = "Inner") {
    Table <- matrix(as.numeric(as.matrix(TableChar)), nrow = dimension[1])
    Table <- fun(Table)
    if (log == TRUE)
        Table <- log(Table)
    else if (log != FALSE)
        Table <- log(Table) / log(log)
    ## If 'fun' reduces the dimensions, then 'Inner' should be
    ## present, and the following two updates modified:
    Inner <- c(Inner, tolower(Inner))
    Inner <- c(Inner, substr(Inner, 1, 1))
    if (length(-which(!is.na(match(names(dimension), Inner)))) > 0) {
        dimension <- dimension[-which(!is.na(match(names(dimension), Inner)))]
        names <- names[-which(!is.na(match(names(names), Inner)))]
        ## 'inversPerm' does not contain 'inner'.
    }
    dimR <- dimension[3:length(dimension)]
    ## Handles 'combinedTreatment', by only if given as 'COLUMS':
    d <- length(inversPerm)
    is.wholenumber <-
        function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if (!all(is.wholenumber(dimR)) & !combinedTreatment) {
        warning(paste("Did you forget 'combinedTreatment = TRUE'",
                      "for 'readAssayTable' ?"), call. = FALSE)
        ## Trying to fix:
        combinedTreatment <- TRUE
    }
    if (combinedTreatment & (2 < d))
        if ((inversPerm[d] == d) & (inversPerm[d-1] == d-1)) {
                if (length(names[d-1][[1]]) == dimension["COLUMN"]) {
                    names[d-1] <- paste(names[d-1], names[d], collapse = ":")
                    dimR[d-1] <- dimension["COLUMN"]
                    dimR[d] <- 1
                }
            }
    if (!all(is.wholenumber(dimR))) {
        warning(paste("Unbalanced design with different number of doses \n",
                      "for different sample can only be defined by columns"),
                call. = FALSE)
    }
    if (all(is.wholenumber(dimR)) & (prod(dim(Table)) == prod(dimR))) {
        dim(Table) <- dimR
        namesR <- NULL
        for (i in 1:length(dimR)) {
            uNames <- unique(names[[i]])
            if (length(uNames) == dimR[i])
                namesR <- append(namesR, list(uNames))
            else
                if (dimR[i] == 1)
                    namesR <- append(namesR, paste(uNames, collapse = " / "))
                else
                    if (dimR[i] * length(uNames) %/% dimR[i] == length(uNames))
                        namesR <- append(
                            namesR,
                            list(apply(
                                matrix(names[[i]],
                                       ncol = length(uNames) %/% dimR[i]),
                                1, FUN = function(x)
                                paste(x, collapse = " / "))))
                    else
                        namesR <- append(namesR, list(NULL))
        }
        dimnames(Table) <- namesR ## lapply(names, FUN = function(x) unique(x))
        Table <- aperm(Table, inversPerm)
    } else {
        warning(paste("Problems with dimensions \n",
                      " - table returned 'as is' in array"), call. = FALSE)
        dimnames(Table) <- dimnames(TableChar)
    }
    ## Dimensions not reduced by combining the two last dimensions!
    return(Table)
}
