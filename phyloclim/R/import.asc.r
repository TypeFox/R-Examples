import.asc <- function (file, type = c("numeric", "factor"), lev = NULL,
                          levnb = 1, labnb = 3)
{
    ## Verifications
    type <- match.arg(type)
    if (substr(file, nchar(file) - 3, nchar(file)) != ".asc")
        stop("not a valid .asc file")
    if ((type != "numeric") & (type != "factor"))
        stop("argument type should be \"factor\" or \"numeric\"")
    if ((type == "numeric") & (!is.null(lev)))
        stop("lev can be specified only when type is \"factor\" ")
    if ((type == "factor") & (length(lev) == 1))
        if (!file.exists(lev))
            stop("lev is not a valid file")

    ## Opens a connection
    zz <- file(file, "r")
    ## reads the header
    nc <- readLines(zz, 1)
    nl <- readLines(zz, 1)
    xll <- readLines(zz, 1)
    yll <- readLines(zz, 1)
    cs <- readLines(zz, 1)
    nas <- readLines(zz, 1)
    cs <- strsplit(cs, " ")
    cs <- as.numeric(cs[[1]][length(cs[[1]])])

    ## The coordinates of the lower left cell:
    ## Are they coordinates of the corner or thos of the center of the cell?
    cornx <- TRUE
    corny <- TRUE

    ## values of xll and yll
    xll <- strsplit(xll, " ")
    if ((xll[[1]][1] == "xllcenter") | (xll[[1]][1] == "XLLCENTER"))
        cornx <- FALSE
    xll <- as.numeric(xll[[1]][length(xll[[1]])])
    yll <- strsplit(yll, " ")
    if ((yll[[1]][1] == "yllcenter") | (xll[[1]][1] == "YLLCENTER"))
        corny <- FALSE
    yll <- as.numeric(yll[[1]][length(yll[[1]])])

    ## code for NAs
    nas <- strsplit(nas, " ")
    nas <- as.numeric(nas[[1]][length(nas[[1]])])

    ## number of columns
    nc <- strsplit(nc, " ")
    nc <- as.numeric(nc[[1]][length(nc[[1]])])

    ## number of rows
    nl <- strsplit(nl, " ")
    nl <- as.numeric(nl[[1]][length(nl[[1]])])



    ## reads the rest of the file
    tmp <- readLines(zz)

    ## and closes the connection
    close(zz)



    ## opens a new connection to a temporary file
    file.create("toto230876.tmp")
    zz <- file("toto230876.tmp", "w")

    ## and write the content of tmp
    writeLines(tmp, zz)

    ## and finally close this connection
    close(zz)

    ## scan this file as usual and remove tmp file
    output <-scan("toto230876.tmp", quiet=TRUE)
    file.remove("toto230876.tmp")

    ## Place the NAs
    output[output == nas] <- NA

    ## and code into a matrix of class "asc"
    output<-matrix(c(as.matrix(output)), ncol=nl)
    output <- output[, ncol(output):1]

    ## In the case of a factor map
    if (type == "factor") {

        ## creates the vector of levels
        ## if is null, created from the matrix
        if (is.null(lev))
            lev <- levels(factor(output))

        ## if not null, should contain the same number of labels
        ## as there are levels
        if (length(lev) > 1) {
            if (length(lev) != length(levels(factor(output))))
                stop("uncorrect length of lev")
        }

        ## if of length one, then read the correspondence
        ## table exported from Arcview
        if (length(lev) == 1) {
            toto <- read.table(lev, header = TRUE, sep = ",")
            toto <- data.frame(lev = toto[, levnb],
                               hihi = rep(1, nrow(toto)),
                               lab = toto[, labnb])
            toto <- toto[order(toto[, 1]), ]
            if (nrow(toto) != nlevels(factor(output)))
                stop("lev is not a valid correspondence table exported from Arcview")
            lev <- as.character(toto[, 3])
        }
        attr(output, "levels") <- lev
    }

    ## rest of the output
    attr(output, "xll") <- xll
    if (cornx)
        attr(output, "xll") <- xll + cs/2
    attr(output, "yll") <- yll
    if (corny)
        attr(output, "yll") <- yll + cs/2
    attr(output, "cellsize") <- cs
    attr(output, "type") <- type
    class(output) <- "asc"
    return(output)
}

