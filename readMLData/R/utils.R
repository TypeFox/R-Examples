split.comma <- function(x)
{
    stopifnot(length(x) == 1)
    strsplit(x, ",")[[1]]
}

join.comma <- function(x)
{
    paste(x, collapse=",")
}

dsSort <- function(dat)
{
    ind <- do.call("order", as.list(dat))
    dat[ind, ]
}

getFields <- function()
{
    colnames(fieldInfo) <- c("field", "type", "usage")
    data.frame(fieldInfo)
}

getPath <- function(dirName)
{
    paste(library(help=readMLData)$path, dirName, sep="/")
}

getIndex <- function(dsList, id)
{
    if (class(id) == "character") {
        ind <- which(dsList$identification == id)
    } else if (is.numeric(id)) {
        ind <- id
    } else {
        stop("identification of a database must be of character or numeric type")
    }
    if (length(ind) == 0) {
        cat("no data set found for", id, "\n")
        return(invisible(NULL))
    }
    if (length(ind) >= 2) {
        cat("several data sets found for", id, "\n")
        return(invisible(NULL))
    }
    ind
}

getType <- function(dat)
{
    type <- rep(NA, times=ncol(dat))
    for (i in 1:length(type))
    {
        col.class <- class(dat[, i])
        stopifnot(length(col.class) == 1)
        col.dat <- dat[, i]
        col.dat <- col.dat[!is.na(col.dat)]
        if (col.class %in% c("character", "factor"))
            type[i] <- as.character(length(unique(col.dat)))
        else if (col.class %in% c("numeric", "double", "integer"))
            type[i] <- "n"
        else
            type[i] <- "o"
    }
    type
}

analyzeData <- function(dat)
{
    n <- ncol(dat)
    col.class <- rep(NA, times=n)
    col.num.unique <- rep(NA, times=n)
    col.min <- rep(NA, times=n)
    col.max <- rep(NA, times=n)
    for (i in seq.int(length.out=ncol(dat))) {
        x <- dat[, i]
        col.class[i] <- class(x)
        col.num.unique[i] <- length(unique(x))
        if (is.numeric(x)) {
            col.min[i] <- min(x)
            col.max[i] <- max(x)
        }
    }
    data.frame(class=col.class, num.unique=col.num.unique, min=col.min, max=col.max)
}

