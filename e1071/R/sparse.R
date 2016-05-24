read.matrix.csr <- function(file, fac = TRUE, ncol = NULL)
{
    l <- strsplit(readLines(file), "[ ]+")

    ## extract y-values, if any
    y <- if (is.na(l[[1]][1]) || length(grep(":",l[[1]][1])))
        NULL
    else
        sapply(l, function(x) x[1])

    ## x-values
    rja <- do.call("rbind",
                   lapply(l, function(x)
                          do.call("rbind",
                                  strsplit(if (is.null(y)) x else x[-1], ":")
                                  )
                          )
                   )
    ja <- as.integer(rja[,1])
    ia <- cumsum(c(1, sapply(l, length) - !is.null(y)))

    max.ja <- max(ja)
    dimension <- c(length(l), if (is.null(ncol)) max.ja else max(ncol, max.ja))
    x = new(getClass("matrix.csr", where = asNamespace("SparseM")),
    ra = as.numeric(rja[,2]), ja = ja,
    ia = as.integer(ia), dimension = as.integer(dimension))
    if (length(y))
        list(x = x, y = if (fac) as.factor(y) else as.numeric(y))
    else x
}

write.matrix.csr <- function (x, file = "out.dat", y = NULL, fac = TRUE) {
    on.exit(sink())

    x <- SparseM::as.matrix.csr(x)
    if (!is.null(y) & (length(y) != nrow(x)))
        stop(paste("Length of y (=", length(y),
                   ") does not match number of rows of x (=",
                   nrow(x), ")!", sep=""))
    sink(file)
    l <- length(x@ra)
    zerocols <- all(x@ja < ncol(x))
    if (!is.null(y) && is.factor(y) && fac)
        y <- as.character(y)
    for (i in 1:nrow(x)) {
        if (!is.null(y)) cat (y[i],"")
        if ((x@ia[i] <= l) && (x@ia[i] < x@ia[i + 1])) {
            for (j in x@ia[i] : (x@ia[i + 1] - 1))
                cat(x@ja[j], ":", x@ra[j], " ", sep="")
            if (zerocols) {
                cat(ncol(x), ":", 0, " ", sep="")
                zerocols <- FALSE
            }
        }
        cat("\n")
    }
}

na.fail.matrix.csr <- function(object, ...) {
    if (any(is.na(object@ra)))
        stop("missing values in object") else return(object)
}






