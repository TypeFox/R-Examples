# --------------------------------------------------------------------
# compare.R
# --------------------------------------------------------------------

.graphicsQCenv <- new.env()

# --------------------------------------------------------------------
#
# compare()
#
# compare will take two qcPlot* objects of the same class and identify
# any differences in the plots they produced.
# compare has the option of deleting the test files and the model files
# once the comparison has been made, or to retain files if they
# are different.
#
# --------------------------------------------------------------------
`compare` <-
function(test,
         control,
         path = test$info$directory,
         erase = c("none", "identical", "files", "all")
         )
{
    # Start warning handler
    assign("graphicsQCWarnings", character(0), envir = .graphicsQCenv)
    # Clear warning handler on exit
    on.exit(rm("graphicsQCWarnings", envir = .graphicsQCenv))

    if (is.character(test) && is.character(control) && test == control) {
        # If auto-detecting in the same dir it's not clear which is test
        # or which is control (or what to do if 3 are found!)
        stop(gettextf("%s and %s paths cannot be the same",
                      sQuote("test"), sQuote("control")),
             domain=NA)
    }
    
    if (!erase[1] %in% c("none", "identical", "files", "all")) {
        warning(gettextf("%s must be one of %s, %s, %s, or %s - %s used",
                         sQuote("erase"), dQuote("none"), 
                         dQuote("files"), dQuote("identical"), ", or ",
                         dQuote("all"), dQuote("none")),
                domain=NA)
        erase <- "none"
    }
    
    test <- getQCResult(test)
    control <- getQCResult(control)

    path <- getValidPath(path)

    if (inherits(test, "list") || inherits(control, "list")) {
        if (inherits(test, "list") && inherits(control, "list")) {
            notYetImplemented("Comparing lists of test and control results")
            RESULT <- mapply(compare, test, list, erase) ## return?
        } else {
            # Note: should also check that the list isn't of length 1
            stop("Cannot have a list compared to one")
        }
    }
    ## first package, then:
    if (inherits(test, "qcPlotFunResult") &&
        inherits(control, "qcPlotFunResult")) {
        return(compareFunOrFile(test, control, path, erase, "Fun"))
    } else if (inherits(test, "qcPlotFileResult") &&
               inherits(control, "qcPlotFileResult")) {
        return(compareFunOrFile(test, control, path, erase, "File"))
    } else if (inherits(test, "qcPlotExprResult") &&
               inherits(control, "qcPlotExprResult")) {
        results <- compareExpr(test, control, path, erase)
        return(results)
    } else {
        # Test and Control are not the same classes!
        stop("test and control are not of the same class")
    }
    results
}

# --------------------------------------------------------------------
#
# compareExpr()
#
# --------------------------------------------------------------------
`compareExpr` <-
function(test, control, path, erase)
{
    filePairs <- getPairs(test, control)
    # names(filepairs) are the filetypes to compare
    results <- lapply(names(filePairs[["test"]]), compareType,
                 filePairs[["control"]], control[["info"]][["directory"]],
                 filePairs[["test"]], test[["info"]][["directory"]],
                 path, erase)
    if (length(filePairs[[1]]) < length(filePairs[[2]])) {
        names(results) <- names(filePairs[[1]])
    } else {
        names(results) <- names(filePairs[[2]])
    }
    results <- compareWarnings(test, control, results)
    results[["unpaired"]] <- filePairs[["unpaired"]]
    info <- list("OS" = .Platform$OS.type, "Rver" =
                 version[["version.string"]], "date" = date(),
                 "call" = paste(deparse(sys.call(sys.parent())),
                   collapse = ""),
                 "path" = normalizePath(path.expand(path)),
                 "testDirectory" = test[["info"]][["directory"]],
                 "controlDirectory" = control[["info"]][["directory"]],
                 "logFilename" = getCompareExprLogFilename(test, control))
    results <- list("info" = info, "testInfo" = test[["info"]],
                    "controlInfo" = control[["info"]], "results" = results)
    class(results) <- "qcCompareExprResult"
    writeXmlCompareExprLog(results)
    results
    # Compare one plotExpr at a time (using natural order of qcresult
    # OR order specified explicitly by log files OR the natural
    # order of the log files from the autodetect).
}

# --------------------------------------------------------------------
#
# compareFunOrFile()
#
# Note: This function isn't to be used by the user. It gets called
#       from compare().
# --------------------------------------------------------------------
`compareFunOrFile` <-
function(test, control, path, erase, type)
{
    # test[[-1]] and control[[-1]] take out info before comparing
    results <- mapply(compareExpr, test[[-1]], control[[-1]], 
                      MoreArgs = list(path = path, erase = erase),
                      SIMPLIFY = FALSE)
    info <- list("OS" = .Platform$OS.type, "Rver" =
                 version[["version.string"]], "date" = date(),
                 "call" = paste(deparse(sys.call(sys.parent())),
                   collapse = ""),
                 "path" = normalizePath(path.expand(path)),
                 "logFilename" =
                 paste(unlist(strsplit(results[[1]][["testInfo"]][[
                       "logFilename"]], "-log.xml")),
                       "-compare", type, "Log.xml", sep = ""),
                 "testLog" = file.path(test[["info"]][["directory"]],
                   test[["info"]][["logFilename"]]),
                 "controlLog" = file.path(control[["info"]][["directory"]],
                   control[["info"]][["logFilename"]]))
    results <- list("info" = info, "results" = results)
    writeXmlCompareTypeLog(results, type)
    class(results) <- paste("qcCompare", type, "Result", sep = "")
    results
}

# --------------------------------------------------------------------
#
# getPairs()
#
# Also sorts unpaired warnings
# --------------------------------------------------------------------
`getPairs` <-
function(test, control)
{
    testFiletypes <- names(test[["plots"]])
    controlFiletypes <- names(control[["plots"]])

    # If the amount of files for a given filetype have different
    # length, put the leftovers in 'unpaired' and cut the group
    # with more files down to size
    # NB: when this is the case, it is likely that one extra plot in
    # the middle of the other plots would cause the rest to fail.

    allFiletypes <- unique(c(testFiletypes, controlFiletypes))
    filetypesToCompare <- controlFiletypes[controlFiletypes %in% testFiletypes]

    controlUnpairedFiletype <- controlFiletypes[!(controlFiletypes %in%
                                                  filetypesToCompare)]
    testUnpairedFiletype <- testFiletypes[!(testFiletypes %in%
                                            filetypesToCompare)]

    # Control filetypes are being compared, everything else is unpaired
    controlPaired <- vector("list", length(controlFiletypes))
    names(controlPaired) <- controlFiletypes
    testPaired <- controlPaired

    testUnpaired <- vector("list", length(testFiletypes))
    names(testUnpaired) <- testFiletypes
    controlUnpaired <- controlPaired
    # First take out all the ones that aren't even paired filetypes
    lapply(controlUnpairedFiletype, function(type) {
        lapply(names(control[["plots"]][[type]]), function (ele)
            if (!is.null(control[["plots"]][[type]][[ele]])) {
                controlUnpaired[[type]][ele] <<- list(control[["plots"]][[
                                                               type]][[ele]])
            })
        if (length(controlUnpaired[[type]]) == 0) {
            controlUnpaired[[type]] <<- NULL
        }
    })
    lapply(testUnpairedFiletype, function(type) {
        lapply(names(test[["plots"]][[type]]), function (ele)
            if (!is.null(test[["plots"]][[type]][[ele]])) {
                testUnpaired[[type]][ele] <<- list(test[["plots"]][[
                                                         type]][[ele]])
            })
        if (length(testUnpaired[[type]]) == 0) {
            testUnpaired[[type]] <<- NULL
        }
    })
    # Take out unpaired files within the paired filetypes
    lapply(filetypesToCompare, function (filetype) {
        testPlotIndices <- seq_along(test[["plots"]][[filetype]][["plot"]])
        controlPlotIndices <- seq_along(control[["plots"]][[
                                        filetype]][["plot"]])
        shortest <- seq_len(min(length(testPlotIndices),
                                length(controlPlotIndices)))
        unpTestPlots <- test[["plots"]][[filetype]][["plot"]][
                                              !(testPlotIndices %in% shortest)]
        unpControlPlots <- control[["plots"]][[filetype]][["plot"]][
                                           !(controlPlotIndices %in% shortest)]
        if (length(unpTestPlots) > 0 || length(unpControlPlots) > 0) {
            warningHandler("length of files to compare are different;",
                           " unpaired files ignored")
        }
        testUnpaired[[filetype]] <<- if (length(unpTestPlots) != 0)
                                         list(plot=unpTestPlots)
        controlUnpaired[[filetype]] <<- if (length(unpControlPlots) != 0)
                                            list(plot=unpControlPlots)
        controlPaired[filetype] <<- list(control[["plots"]][[filetype]][[
                                                  "plot"]][shortest])
        testPaired[filetype] <<- list(test[["plots"]][[filetype]][[
                                           "plot"]][shortest])
        })
    # If the unpaireds are just list() or blank, change them to NULL
    # for consistency when reading the XML file later on
    if (length(testUnpaired) == 0) testUnpaired <- NULL
    if (length(controlUnpaired) == 0) controlUnpaired <- NULL
    return(list("test" = testPaired, "control" = controlPaired, "unpaired" =
         list("test"=testUnpaired, "control"=controlUnpaired)))
}

# --------------------------------------------------------------------
#
# compareWarnings()
#
# Only compares the warnings/errors in the filetypes being compared. Ones
# in completely unpaired filetypes are dealt with by getPairs().
#
# --------------------------------------------------------------------
`compareWarnings` <-
function(test, control, results)
{
    types <- lapply(results, length)
    # The types of length 0 are completely unpaired so we don't compare them
    # here.
    for (type in names(which(types > 0))) {
        contWarns <- control[["plots"]][[type]][["warnings"]]
        testWarns <- test[["plots"]][[type]][["warnings"]]
        if ((!is.null(contWarns) || !is.null(testWarns)) &&
            length(contWarns) != length(testWarns) ||
            any(contWarns != testWarns)) {
            results[[type]][["controlWarnings"]] <- contWarns
            results[[type]][["testWarnings"]] <- testWarns
        }
        contError <- control[["plots"]][[type]][["error"]]
        testError <- test[["plots"]][[type]][["error"]]
        if ((!is.null(contError) || !is.null(testError)) &&
            length(contError) != length(testError) ||
            any(contError != testError)) {
            results[[type]][["controlError"]] <- contError
            results[[type]][["testError"]] <- testError
        }
    }
    results
}

## Note: To extend this to allow for more filetypes, just add a function
## "compareNEWTYPE". Once it has been added to valid types in makeplots, that
## should be all the changes required here.
# --------------------------------------------------------------------
#
# compareType()
#
# --------------------------------------------------------------------
`compareType` <-
function(filetype, control, controlPath, test, testPath, path, erase)
{
    # pastes "compare" and 'filetype' to call the appropriate function;
    # pastes 'path' and 'filename' for control and test groups respectively;
    # passes IM capability if it's supported and a diff plot is required
    if (length(test[[filetype]]) == 0) {
        return(NULL)
    } else {
        result <- mapply(paste("compare", toupper(filetype), sep = ""),
                         file.path(testPath, test[[filetype]]),
                         file.path(controlPath, control[[filetype]]),
                         hasIM() && filetype %in% getSupportedIMFormats() &&
                         any(erase == c("none", "identical")),
                         path,
                         SIMPLIFY = FALSE)
        names(result) <- NULL
        return(result)
    }
    ## now remove files?
}

# --------------------------------------------------------------------
#
# comparePDF()
#
# --------------------------------------------------------------------
`comparePDF` <-
function(file1, file2, useIM, diffPlotPath)
{
    # For PDF just compare and ignore the first 6 lines (creationdate/moddate)
    diffName <- getDiffName(file1, file2)
    diffFileName <- paste(diffName, ".diff", sep = "")
    diffPlotName <- paste(diffName, ".png", sep = "")
    diffFilePath <- file.path(diffPlotPath, diffFileName)
    diffResult <- GNUdiff(file1, file2, diffFilePath)
    if (diffResult == "different" && length(readLines(diffFilePath, n = 7))
                                                                        > 6) {
        # There is a true difference, not just the dates/times
        if (useIM) {
            diffPlot <- file.path(diffPlotPath, diffPlotName)
            makeIMDiffPlot(file1, file2, diffPlot)
        }
    } else {
        # Files are the same or just the dates/times were different
        file.remove(diffFilePath)
        diffFileName <- diffPlotName <- NULL
        diffResult <- "identical"
    }
    return(list(controlFile=file2, testFile=file1, result=diffResult,
                diffFile=diffFileName, diffPlot=diffPlotName))
}

# --------------------------------------------------------------------
#
# comparePS()
#
# --------------------------------------------------------------------
`comparePS` <-
function(file1, file2, useIM, diffPlotPath)
{
    diffName <- getDiffName(file1, file2)
    diffFileName <- paste(diffName, ".diff", sep = "")
    diffPlotName <- paste(diffName, ".png", sep = "")
    diffFilePath <- file.path(diffPlotPath, diffFileName)
    diffResult <- GNUdiff(file1, file2, diffFilePath)
    if (diffResult == "different") {
        if (useIM) {
            diffPlot <- file.path(diffPlotPath, diffPlotName)
            makeIMDiffPlot(file1, file2, diffPlot)
        }
    } else {
        diffFileName <- diffPlotName <- NULL
    }
    return(list(controlFile=file2, testFile=file1, result=diffResult,
                                diffFile=diffFileName, diffPlot=diffPlotName))
}

# --------------------------------------------------------------------
#
# comparePNG()
#
# --------------------------------------------------------------------
`comparePNG` <-
function(file1, file2, useIM, diffPlotPath)
{
    diffName <- getDiffName(file1, file2)
    diffPlotName <- NULL
    diffResult <- GNUdiff(file1, file2)
    if (useIM && diffResult == "different") {
        diffPlotName <- paste(diffName, ".png", sep = "")
        diffPlot <- file.path(diffPlotPath, diffPlotName)
        makeIMDiffPlot(file1, file2, diffPlot)
    }
    return(list(controlFile=file2, testFile=file1, result=diffResult,
                                diffFile=NULL, diffPlot=diffPlotName))
}

# --------------------------------------------------------------------
#
# compareBMP()
#
# --------------------------------------------------------------------
`compareBMP` <-
function(file1, file2, useIM, diffPlotPath)
{
    diffName <- getDiffName(file1, file2)
    diffPlotName <- paste(diffName, ".png", sep = "")
    diffResult <- GNUdiff(file1, file2)
    if (useIM && diffResult == "different") {
        diffPlot <- file.path(diffPlotPath, diffPlotName)
        makeIMDiffPlot(file1, file2, diffPlot)
    }
    return(list(controlFile=file2, testFile=file1, result=diffResult,
                                diffFile=NULL, diffPlot=diffPlotName))
}

# --------------------------------------------------------------------
#
# getDiffName()
#
# --------------------------------------------------------------------
`getDiffName` <-
function(file1, file2)
{
    set1 <- basename(file1)
    set2 <- basename(file2)
    paste(gsub("[.]", "-", set1), "+",  gsub("[.]", "-", set2), sep = "")
}

# --------------------------------------------------------------------
#
# GNUdiff()
#
# --------------------------------------------------------------------
`GNUdiff` <-
function(file1, file2, outDiffFile = NULL)
{
    #diffArgs = "-q", intern = FALSE)
    ## *nix only? system() + exit status
    ## This requires a bit more work for Windows support
    # However, XML does not appear to be supported on Mac,
    # and dodgy on Windows (and Sxslt).
    redirectOutput <- ""
    if (!is.null(outDiffFile)) {
        redirectOutput = paste(">", outDiffFile)
    }
    diffResult <- gqcShell(paste("diff", file1, file2, redirectOutput),
                           ignore.stderr=TRUE)
    if (diffResult == 0) {
        # Delete empty diff file
        if (!is.null(outDiffFile)) {
            file.remove(outDiffFile)
        }
        return("identical")
    } else {
        # If one of the files doesn't exist, the .diff file will be empty
        if (!file.exists(file1)) {
            warning(gettextf("file %s not found; marked as different",
                             file1),
                    domain=NA)
            if (!is.null(outDiffFile)) {
                file.remove(outDiffFile)
            }
        }
        if (!file.exists(file2)) {
            warning(gettextf("file %s not found; marked as different",
                             file2),
                    domain=NA)
            if (!is.null(outDiffFile)) {
                file.remove(outDiffFile)
            }
        }
        return("different")
    }

}

# --------------------------------------------------------------------
#
# makeIMDiffPlot()
#
# --------------------------------------------------------------------
`makeIMDiffPlot` <-
function(file1, file2, newFilename)
{
    gqcShell(paste("compare", file1, file2, newFilename))
}

# --------------------------------------------------------------------
#
# hasDiff()
#
# --------------------------------------------------------------------
`hasDiff` <-
function()
{
    length(grep("GNU diffutils", try(gqcShell("diff -v",
                                              intern = TRUE)[1]))) > 0
}

# --------------------------------------------------------------------
#
# hasIM()
#
# --------------------------------------------------------------------
`hasIM` <-
function()
{
    length(grep("ImageMagick", try(gqcShell("compare -version",
                                            intern = TRUE)[1]))) > 0
}

# --------------------------------------------------------------------
#
# getSupportedIMFormats()
#
# --------------------------------------------------------------------
`getSupportedIMFormats` <-
function()
{
    supportedFormats <- character(0)
    formats <- gqcShell("identify -list Format", intern = TRUE,
                        ignore.stderr = TRUE) ## This command may have trouble.
    bmpLine <- grep("Microsoft Windows bitmap image$", formats, value = TRUE,
                    perl = TRUE, useBytes = TRUE)
    pdfLine <- grep("Portable Document Format$", formats, value = TRUE,
                    perl = TRUE, useBytes = TRUE)
    pngLine <- grep("Portable Network Graphics", formats, value = TRUE,
                    perl = TRUE, useBytes = TRUE)
    psLine <- grep("PostScript$", formats, value = TRUE, perl = TRUE,
                   useBytes = TRUE)
    if (length(bmpLine > 0) && grep("[ ]r[^ ][^ ][ ]", bmpLine) > 0) {
        supportedFormats <- c(supportedFormats, "bmp")
    }
    if (length(pdfLine > 0) && grep("[ ]r[^ ][^ ][ ]", pdfLine) > 0) {
        supportedFormats <- c(supportedFormats, "pdf")
    }
    if (length(pngLine > 0) && grep("[ ]rw[^ ][ ]", pngLine) > 0) {
        supportedFormats <- c(supportedFormats, "png")
    }
    if (length(psLine > 0) && grep("[ ]r[^ ][^ ][ ]", psLine) > 0) {
        supportedFormats <- c(supportedFormats, "ps")
    }
    return(supportedFormats)
}

# --------------------------------------------------------------------
#
# mergeList()
#
# Merges elements in named lists which have repeated names
# --------------------------------------------------------------------
`mergeList` <-
function(x)
{
    if (any(duplicated(names(x)))) {
        tags <- unique(names(x))
        output <- lapply(tags, function (tag)
                                   as.vector(unlist(x[names(x) == tag])))
        names(output) <- tags
        return(output)
    } else {
        return(x)
    }
}

# --------------------------------------------------------------------
#
# warningHandler()
#
# --------------------------------------------------------------------
`warningHandler` <-
function(...)
{
    stringWarning <- paste(..., sep = "")
    # Only show warnings we haven't seen before
    qcWarnings <- get("graphicsQCWarnings", .graphicsQCenv)
    if (!stringWarning %in% qcWarnings) {
        assign("graphicsQCWarnings", c(qcWarnings, stringWarning),
               envir = .graphicsQCenv)
        warning(stringWarning, call. = FALSE)
    }
}

# --------------------------------------------------------------------
#
# print.qcCompareExprResult()
#
# --------------------------------------------------------------------
`print.qcCompareExprResult` <-
function (x, ...)
{
    cat("qcCompareExpr Result:\n")
    cat("Call:\n", x[["info"]][["call"]], "\n")
    firstColumn <- format(c("", "R version: ", "Directory: ", "Filename:  "))
    testColumn <- c("             Test", x[["testInfo"]][["Rver"]],
                    shortenPath(x[["testInfo"]][["directory"]]),
                    x[["testInfo"]][["logFilename"]])
    controlColumn <- c("           Control",
                       x[["controlInfo"]][["Rver"]],
                       shortenPath(x[["controlInfo"]][["directory"]]),
                       x[["controlInfo"]][["logFilename"]])
    resultColumn <- c(" Results", rep("", 3))
    lengths <- sapply(x[["results"]],
        function(x) {
            if (is.null(names(x))) {
                length(x)
            } else {
                length(x[names(x) == ""])
            }
        })
    lengths <- lengths[names(lengths) != "unpaired"]
    firstCol <- c("Format:", unlist(lapply(names(lengths), function(type)
                              c(type, rep("", if (lengths[type] - 1 > 0)
                                                  lengths[type] - 1 else 0)))))
    testCol <- controlCol <- resultCol <- character(sum(lengths) + 1)
    testCol[1] <- controlCol[1] <- resultCol[1] <- ""
    i = 2
    for (type in names(lengths)) {
        for (j in seq_len(lengths[type])) {
            testCol[i] <- shortenPath(x[["results"]][[type]][[j]][[
                                           "testFile"]])
            controlCol[i] <- shortenPath(x[["results"]][[type]][[j]][[
                                              "controlFile"]])
            resultCol[i] <- shortenPath(x[["results"]][[type]][[j]][[
                                              "result"]])
            i <- i + 1
        }
        if (lengths[type] == 0) {
            testCol[i] <- controlCol[i] <- "none"
            resultCol[i] <- ""
            i <- i + 1
        }
    }
    mapply(
        function(firstColumn, testColumn, controlColumn, resultColumn) {
            cat(firstColumn, testColumn, controlColumn, resultColumn, "\n",
                sep = "")
        }, format(c(firstColumn, firstCol), width = 11),
           format(c(testColumn, testCol), width = 30),
           format(c(controlColumn, controlCol), width = 30),
           format(c(resultColumn, resultCol), width = 9))
}

# --------------------------------------------------------------------
#
# shortenPath()
#
# Used by print.qcCompareExprResult()
# --------------------------------------------------------------------
`shortenPath` <-
function(path)
{
    lengthPath <- nchar(path)
    if (lengthPath > 29) {
        path <- unlist(strsplit(path, ""))[(lengthPath - 25):lengthPath]
        path <- paste("...", paste(path, collapse = ""), sep = "")
    }
    path
}

# --------------------------------------------------------------------
gqcShell <- function(cmd, intern=FALSE, ignore.stderr=FALSE) {
    if (.Platform$OS.type == "windows") {
        diffResult <- shell(cmd, intern=intern, ignore.stderr = ignore.stderr)
    } else {
        diffResult <- system(cmd, intern=intern, ignore.stderr = ignore.stderr)
    }    
}

