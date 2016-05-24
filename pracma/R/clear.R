##
##  c l e a r . R
##


clear <- function(lst) {
    if (missing(lst))
        lst <- ls(name = .GlobalEnv)
    if (!is.character(lst))
        stop("Argument must be empty or a character vector.")
    rm(list = lst, envir = globalenv())
    # capture.output(gc())
    null <- gc()
}


who <- function() ls(name = .GlobalEnv)


whos <- function() {
    envir <- parent.frame()
    lslist <- ls(envir)
    if (isempty(lslist))
        return(invisible(NULL))

    m <- max(nchar(lslist))
    for (item in lslist) {
        itemObj   <- eval(parse(text = item), parent.frame())
        itemClass <- class(itemObj)
        itemSize  <- object.size(itemObj)
        itemDim   <- paste(dim(itemObj), collapse="x")
        if (itemDim == '') itemDim <- length(itemObj)

        itemSize <- as.numeric(itemSize)
        if (itemSize < 1024) itemSize <- paste(itemSize, "Byte")
        else if (itemSize >= 1024 & itemSize < 1024*1024)
            itemSize <- paste(round(itemSize/1024, 1), "KB")
        else
            itemSize <- paste(round(itemSize/1024/1024, 1), "MB")

        format(cat( item, blanks(m - nchar(item) + 2),
                    itemClass, ", ",
                    itemDim, ", ",
                    itemSize, "\n", 
                    sep=""), justify="centre")
    }
    cat("\n")
    invisible(lslist)
}


what <- function(dname = getwd()) {  # , fexp = "*.R"
    if (is.na(file.info(dname)$isdir)) {
        cat("Argument '", dname, "' is not a known directory.\n", sep = '')
    } else if (file.info(dname)$isdir) {
        fnames <- list.files(dname)
        if (isempty(fnames)) {
            cat("No files in Directory ", dname, ".\n\n", sep = '')
        } else {
            cat("Files in Directory ", dname, ":\n\n", sep = '')
            for (fname in fnames) {
                gname <- paste(dname, fname, sep = "/")
                if (!file.info(gname)$isdir) {
                    cat(fname, "\n")
                }
            }
        }
    } else {
        cat("Argument '", dname, "' is not a directory.\n", sep = '')
    }
    invisible(NULL)
}


cd <- function(dname) {
    if (missing(dname)) {
        dname <- getwd()
    } else {
        setwd(dname)
        dname <- getwd()
    }
    return(dname)
}


pwd <- function()
    getwd()

