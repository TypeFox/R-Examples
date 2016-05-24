# --------------------------------------------------------------------
# makeplots.R
# --------------------------------------------------------------------

validFiletypes <- function() {
    filetypes <- c("pdf", "png", "ps", "bmp")
    if (.Platform$OS.type != "windows") {
        filetypes <- filetypes[-4]
    }
    filetypes
}

# --------------------------------------------------------------------
#
# plotExpr()
#
# plotExpr will take an R expression (or expressions) and produce
# plots in specific file formats of the intended plot.
#
# --------------------------------------------------------------------
`plotExpr` <-
function(expr, # character vector
               # R expression(s)
               # May be a list BUT if it is
               # we just flatten it to a vector.
         filetype = NULL, # character vector
                          # (valid) file formats
         path = NULL, # char length 1
                      # directory to create files in
         prefix = "graphicsqc", # char length 1
                                # file prefix
         clear = FALSE # boolean, clear any files we might make
        )
{
    # These warnings not in evalPlotCode are not recorded in the XML.
    if (is.list(expr))
        expr <- unlist(expr)

    # Testing filetype is valid
    if (is.null(filetype)) {
        filetype <- validFiletypes()
    } else {
        filetype <- getValidFiletypes(filetype)
    }
    fileExtension <- paste(".", filetype, sep = "")

    # Testing file prefix is valid
    # Prefix must be length 1
    if (length(prefix) == 1) {
        prefix <- as.character(prefix)
    } else {
        warning(gettextf("%s should be a character vector of length 1, only the first used (%s)",
                         sQuote("prefix"), prefix[1]),
                domain=NA)
        prefix <- prefix[1]
    }
    wd <- getwd()
    on.exit(setwd(wd))
    # Get valid path
    path <- getValidPath(path)

    # Check we are not going to overwrite any files
    filenamePattern <- paste("^(", prefix, "-", "[0-9]+[.](",
                   paste(filetype, collapse = "|"),
                   ")|", prefix, "-log[.]xml)$", sep = "")
                   # ie "^(prefix-[0-9]+[.](pdf|ps)|prefix-log[.]xml)$"
                   # Will match (prefix = "prefix"; filetype = c("pdf", "ps")):
                   # prefix-123.pdf, prefix-log.xml, prefix-1.ps
                   # Will NOT match:
                   # prefix-log.pdf, prefix-123.xml, prefix-.ps
    currentFilenames <- list.files(path, filenamePattern)
    # 'Clear' files we might make if we are told to
    if (length(clear) > 1) {
        warning(gettextf("%s has more than one element: only the first used",
                         sQuote("clear")),
                domain=NA)
        clear <- clear[1]
    }
    if (is.logical(clear) && !is.na(clear)) {
        if (clear && length(currentFilenames) > 0) {
            if (any(!file.remove(file.path(path, currentFilenames)))) {
                stop("tried to clear but failed to delete files")
            }
        } else if (!clear && length(currentFilenames) > 0) {
            stop(gettextf("files of intended filename already exist in %s",
                          sQuote("path")),
                 call. = FALSE, domain=NA)
        }
    } else {
        stop(gettextf("%s must be either TRUE or FALSE", sQuote("clear")),
             domain=NA)
    }

    filenameFormat <- paste(prefix, "-%d", sep = "")
    setwd(path)
    evalResult <- lapply(filetype, evalPlotCode, expr, filenameFormat)
    names(evalResult) <- filetype
    # ---------- Get info of results ----------
    # (We only created files if none already existed)
    blankImageSizes <- getBlankImageSizes()
    if (any(is.na(blankImageSizes))) {
        generateBlankImages()
        blankImageSizes <- getBlankImageSizes()
        if (any(is.na(blankImageSizes))) {
            warning("blank image details could not be obtained, so blank files will not be removed")
        }
    }

    filenames <- list.files(getwd(), filenamePattern)
    # Remove blanks
    wereRemoved <- sapply(filenames, removeIfBlank, blankImageSizes,
                          USE.NAMES = FALSE)
    if (any(!wereRemoved)) {
           warning("some blank images could not be removed")
    }
    filenames <- list.files(getwd(), filenamePattern) ##full.names = TRUE
    plots <- lapply(filetype,
                function(filetype) {
                    plots <- grep(filetype, filenames, value = TRUE)
                    if (length(plots) > 0) plots else NULL
                })
    names(plots) <- filetype
    lapply(filetype, function(type)
                         evalResult[[type]]["plot"] <<- list(plots[[type]]))

    call <- deparse(sys.call(sys.parent()))
    if (length(grep("^plotExpr", call)) == 0) {
        call <- deparse(sys.call(sys.parent(2)))
        if (length(grep("^plotFile", call)) > 0) {
            call <- "Called from plotFile()"
        }
        if (length(grep("^plotFunction", call)) > 0) {
            call <- "Called from plotFunction()"
        }
    }
    info <- list("OS" = .Platform$OS.type, "Rver" =
                 version[["version.string"]], "date" = Sys.time(),
                 "call" = paste(call, collapse = "\n"),
                 ## at some point deparse(width.cutoff) might need to be raised
                 "directory" = getwd(), "logFilename" =
                 paste(prefix, "-log.xml", sep = ""))
                 # using getwd() here forces expansion
                 # of directory (ie, expands "./" or "~/")
    results <- list("info" = info, "plots" = evalResult)

    writeXmlPlotExprLog(results)
    class(results) <- "qcPlotExprResult"
    # ---------- return results ----------
    invisible(results)
}

# --------------------------------------------------------------------
#
# getValidPath()
#
# --------------------------------------------------------------------
`getValidPath` <-
function(path)
{
    if (length(path) == 0) {
        warning("no path given: the path has been set to your current working directory", call. = FALSE)
        path <- getwd()
    } else if (length(path) == 1) {
        isDir <- file.info(path)$isdir
        isCreated <- TRUE
        if (is.na(isDir)) {
            isCreated <- dir.create(path, showWarnings = FALSE)
        } else if (!isDir || !isCreated) {
            stop(gettextf("could not create directory %s", dQuote(path)),
                 call. = FALSE, domain=NA)
        }
    } else {
        warning(gettextf("the given path has more than one element: only the first used in %s",
                         dQuote(path[1])),
                call. = FALSE, domain=NA)
        path <- getValidPath(path[1])
    }
    path
}

# --------------------------------------------------------------------
#
# evalPlotCode()
#
# --------------------------------------------------------------------
`evalPlotCode` <-
function(filetype, expr, filenameFormat)
{
    if (filetype == "ps") {
        fileExtension <- ".ps"
        filetype <- "postscript"
    } else {
        fileExtension <- paste(".", filetype, sep = "")
    }
    if (any(filetype == c("pdf", "postscript"))) {
        do.call(filetype, list(paste(filenameFormat, fileExtension,
                                     sep = ""), onefile = FALSE))
    } else {
        do.call(filetype, list(paste(filenameFormat, fileExtension,
                                     sep = "")))
    }
    warns <- NULL
    # Reset last error message to ""
    tryCatch(stop(), error = function(e) {})
    tryCatch(withCallingHandlers(eval(if (is.language(expr)) expr else
                                                         parse(text = expr)),
                    warning = function(w) {
                      warns <<- c(warns, paste(conditionMessage(w)));
                      invokeRestart("muffleWarning")
                    }), error = function(e) {})
    error <- geterrmessage() # There can only be one error as we stop
                             # evaluating when we hit an error
    # Strip any trailing newlines as XML won't read them in properly later
    # (will write them, but won't read the last \n's (as normal or CData))
    warns <- sub("[\n]*$", "", warns)
    error <- sub("[\n]*$", "", error)
    if (length(error) == 0 || error == "") {
        error <- NULL
    }
    if (length(warns) == 0) {
        warns <- NULL
    }
    dev.off()
    return(list("warnings" = warns, "error" = error))
}

# --------------------------------------------------------------------
#
# getValidFiletypes()
#
# --------------------------------------------------------------------
`getValidFiletypes` <-
function(filetypes)
{
    filetypes <- tolower(filetypes)
    validFiletypes <- validFiletypes()

    # check for duplications
    if (any(duplicated(filetypes))) {
        warning(gettextf("duplicated filetypes: %s duplication ignored",
                         paste(dQuote(filetypes[duplicated(filetypes)]),
                               collapse = ", ")),
                call. = FALSE, domain=NA)
        filetypes <- filetypes[!duplicated(filetypes)]
    }

    # check given filetypes against valid filetypes
    invalidTypes <- !filetypes %in% validFiletypes
    if (any(invalidTypes)) {
       if (any(filetypes[invalidTypes] %in% "bmp")) {
           warning("sorry, BMP format only supported in Windows",
                   call. = FALSE)
       }
       warning(gettextf("invalid filetype(s) given: %s ignored",
                        paste(dQuote(filetypes[invalidTypes]),
                              collapse = ", ")),
               call. = FALSE, domain=NA)
    }

    if (length(filetypes[!invalidTypes]) > 0) {
        return(filetypes[!invalidTypes])
    } else {
        stop("no valid filetypes given", call. = FALSE)
    }
}

# --------------------------------------------------------------------
#
# plotFile()
#
# --------------------------------------------------------------------
`plotFile` <-
function(filename, # character vector
                   # R expression(s)
         filetype = NULL, # character vector
                          # (valid) file formats
         path = NULL, # char length 1
                      # directory to create files in
         prefix = basename(filename), # char length 1
                                      # file prefix
         clear = FALSE
         )
{
    ## Test if files exist first?
    if (length(filename) != length(prefix)) {
        stop(gettextf("%s must be the same length as %s",
                      sQuote(prefix), sQuote(filename)),
             domain=NA)
    }
    if (length(grep(.Platform$file.sep, prefix)) > 0) {
        stop(gettextf("%s cannot contain the system file separator",
                      sQuote(prefix)),
             domain=NA)
    }
    path <- getValidPath(path)
    expr <- lapply(filename, readLines)
    fileMapplyResult <- mapply(plotExpr, expr = expr, prefix = prefix,
                              MoreArgs = list(filetype = filetype, path = path,
                              clear = clear), SIMPLIFY = FALSE)
    names(fileMapplyResult) <- NULL
    if (length(prefix) > 1) {
        # No warning - there is a note in the help file
        filePrefix <- prefix[1]
    } else {
        filePrefix <- prefix
    }
    
    path <- normalizePath(path.expand(path))
    info <- list("OS" = .Platform$OS.type, "Rver" =
                 version[["version.string"]], "date" = Sys.time(),
                 "call" = paste(deparse(sys.call(sys.parent())),
                   collapse = "\n"),
                 "directory" = path,
                 "logFilename" = paste(filePrefix, "-fileLog.xml", sep = ""))
    
    # Note: plotExpr XML file gets written in the call to plotExpr
    writeXmlPlotTypeLog(prefix, info, "file")
    fileResult <- list("info" = info, "results" = fileMapplyResult)
    class(fileResult) <- c("qcPlotFileResult")
    fileResult
}

# --------------------------------------------------------------------
#
# plotFunction()
# ## setRNG?
# --------------------------------------------------------------------
`plotFunction` <-
function(fun, # character vector
              # R expression(s)
         filetype = NULL, # character vector
                          # (valid) file formats
         path = NULL, # char length 1
                      # directory to create files in
         prefix = fun, # char length 1
                       # file prefix
         clear = FALSE
         )
{
    ## if not paste, then something like...
    # funs<-lapply(fun, function(x) {
       #                 substitute(do.call(example,list(x))
       # (gets unusual behaviour without paste..?)
    ## can also check if it's a function? .. is.function(match.fun(..
    if (!is.character(fun)) {
        fun <- deparse(substitute(fun))
        prefix <- fun
    }
    path <- getValidPath(path)
    # Check that the prefixes are unique (otherwise we will die half-way
    # through because we will try to overwrite)
    if (any(duplicated(prefix))) {
        stop(gettextf("all values for %s must be unique", sQuote("prefix")),
             domain=NA)
    }

    if (length(fun) != length(prefix)) {
        stop(gettextf("%s must be the same length as %s",
                      sQuote("prefix"), sQuote("fun")),
             domain=NA)
    }
 
    funs <- paste("example(", fun, ", echo = FALSE, setRNG = TRUE)", sep = "")
    funMapplyResult <- mapply(plotExpr, expr = funs, prefix = prefix,
                              MoreArgs = list(filetype = filetype, path = path,
                              clear = clear), SIMPLIFY = FALSE)
    names(funMapplyResult) <- NULL
    if (length(prefix) > 1) {
        # No warning - there is a note in the help file
        filePrefix <- prefix[1]
    } else {
        filePrefix <- prefix
    }
    
    path <- normalizePath(path.expand(path))
    info <- list("OS" = .Platform$OS.type, "Rver" =
                 version[["version.string"]], "date" = Sys.time(),
                 "call" = paste(deparse(sys.call(sys.parent())),
                   collapse = "\n"),
                 "directory" = path,
                 "logFilename" = paste(filePrefix, "-funLog.xml", sep = ""))
    
    writeXmlPlotTypeLog(prefix, info, "fun")
    funResult <- list("info" = info, "results" = funMapplyResult)
    class(funResult) <- c("qcPlotFunResult")
    funResult
}

# --------------------------------------------------------------------
#
# plotPackage()
#
# --------------------------------------------------------------------
`plotPackage` <-
function(package)
{
    ## take 'best effort' approach? try get tests/demo if can and use them..
    if (!do.call(require, list(package))) {
        warning(gettextf("failed to load package %s", dQuote(package)),
                domain=NA)
    } # now package is loaded
    notYetImplemented("producing plots from the name of a package")
}

# --------------------------------------------------------------------
#
# generateBlankImages()
#
# --------------------------------------------------------------------
`generateBlankImages` <-
function()
{
    tempDir <- tempdir()
    # Only pdf and ps blank images are made as the filesize for bmp and png
    # blanks are 0
    pdf(file.path(tempDir, "blankPDF.pdf"),
        onefile = FALSE)
    dev.off()
    postscript(file.path(tempDir, "blankPS.ps"),
               onefile = FALSE)
    dev.off()
    invisible()
}

# --------------------------------------------------------------------
#
# getBlankImageSizes()
#
# --------------------------------------------------------------------
`getBlankImageSizes` <-
function()
{
    sizes <- file.info(file.path(tempdir(),
                                 c("blankPDF.pdf", "blankPS.ps")))[,1]
    names(sizes) <- c("pdf", "ps")
    sizes
}

# --------------------------------------------------------------------
#
# removeIfBlank()
#
# --------------------------------------------------------------------
`removeIfBlank` <-
function(filename, blankImageSizes)
{
    filesize <- file.info(filename)[,1]
    if (filesize == 0 || length(grep(".*[.]pdf$", filename)) > 0 &&
        filesize == blankImageSizes["pdf"] ||
        length(grep(".*[.]ps$", filename)) > 0 && filesize ==
                                               blankImageSizes["ps"]) {
        if (!file.remove(filename)) {
            return(FALSE)
        }
    }
    return(TRUE)
}

# --------------------------------------------------------------------
#
# print.qcPlotExprResult()
#
# --------------------------------------------------------------------
`print.qcPlotExprResult` <-
function (x, ...)
{
    cat("plotExpr Result:\n")
    cat("Call:      ",
        gsub("\n", "\n           ", x[["info"]][["call"]]), "\n",
        "R version: ", x[["info"]][["Rver"]], "\n",
        "Directory: ", x[["info"]][["directory"]], "\n",
        "Filename:  ", x[["info"]][["logFilename"]], "\n", sep="")
    cat("Formats:\n")
    for (format in names(x[["plots"]])) {
        if (is.null(x[["plots"]][[format]][["plot"]])) {
            plots = "none"
        } else {
            plots = x[["plots"]][[format]][["plot"]]
        }
        cat(" ", format, ": Plots: ")
        cat(plots, sep = ", ")
        cat("\n")

        if (!is.null(x[["plots"]][[format]][["warnings"]])) {
            cat("    Warnings: ")
            cat(x[["plots"]][[format]][["warnings"]],
                sep = "\n              ")
            cat("\n")
        }
        if (!is.null(x[["plots"]][[format]][["error"]])) {
            cat("    Error: ", x[["plots"]][[format]][["error"]], "\n")
        }
    }
}


# --------------------------------------------------------------------
#
# notYetImplemented()
#
# --------------------------------------------------------------------
`notYetImplemented` <-
function(feature)
{
    stop(gettextf("Sorry, %s is not yet implemented", feature),
         domain=NA)
}
