# --------------------------------------------------------------------
# report.R
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#
# writeReport()
#
# writeReport will produce a HTML output of any differences found
# from the compare function.
#
# 'comparison' is either a character vector specifying the path to the log
# file to report on, or the folder where it will report on all the log files,
# or an R object to create the report.
# Reports are made for all log files that the given log file might refer to.
# A character vector is returned containing paths to the log files (with only
# the highest classed ones returned when a folder is given).
#
# Note: Currently, if the R object is given, the log file and sub-class
# log files must exist (xsltApplyStyleSheet requires a filename or a string
# containing the doc)
# --------------------------------------------------------------------
`writeReport` <-
function(qcResult, browse = TRUE, xslStyleSheets = NULL)
{
    #SxsltInitializationFunction()
    ## Note: There is a call to browser() in addXSLTFunctions...
    # -- can't just give logToHTML even if it is defined outside - have
    # to give the function definition here?
    Sxslt::addXSLTFunctions(logToHTML = function(...) logNameToHTML(file.path(...)),
                     getCompareExprName =
                       function(logWithPath) {
                           strsplit(basename(logWithPath),
                                    "-compareExprLog.xml")[[1]]
                       }
                    )
    xslStyles <- list("plotExprStyleSheet" =
                      system.file("xsl", "plotExpr.xsl",
                                  package = "graphicsQC"),
                      "plotFunAndFileStyleSheet" =
                      system.file("xsl", "plotFunAndFile.xsl",
                                  package = "graphicsQC"),
                      "compareExprStyleSheet" =
                      system.file("xsl", "compareExpr.xsl",
                                  package = "graphicsQC"),
                      "compareFunAndFileStyleSheet" =
                      system.file("xsl", "compareFunAndFile.xsl",
                                  package = "graphicsQC"))
    xslStyles[names(xslStyleSheets)] <- xslStyleSheets

    if (is.character(qcResult)) {
        # Either a folder or a file.
        fileInfo <- file.info(qcResult)
        if (is.na(fileInfo$isdir)) {
            stop(gettextf("file %s not found - no log file for this created",
                          dQuote(qcResult)),
                 call. = FALSE, domain=NA)
        } else {
            # getQCResult also does auto-detect
            qcResult <- getQCResult(qcResult, "all")
        }
    }

    # Now we have the highest level R object we want to report on.
    filename <- report(qcResult, xslStyles)
    if (browse) {
        browseURL(filename)
    } 
    filename
}

# --------------------------------------------------------------------
#
# report()
#
# --------------------------------------------------------------------
`report` <-
function(qcResult, xslStyles)
{
    UseMethod("report")
}

# --------------------------------------------------------------------
#
# report.qcPlotExprResult()
#
# --------------------------------------------------------------------
`report.qcPlotExprResult` <-
function(qcResult, xslStyles)
{
    plotExprPath <- file.path(qcResult[["info"]][["directory"]],
                              qcResult[["info"]][["logFilename"]])
    plotExpr <- Sxslt::xsltApplyStyleSheet(plotExprPath,
                                    xslStyles[["plotExprStyleSheet"]])
    logName <- logNameToHTML(plotExprPath)
    saveXML(plotExpr$doc, file = logName)
    logName
}

# --------------------------------------------------------------------
#
# report.qcPlotFileResult() and report.qcPlotFunResult()
#
# --------------------------------------------------------------------
`report.qcPlotFileResult` <- `report.qcPlotFunResult` <- 
function(qcResult, xslStyles)
{
    plotFPath <- file.path(qcResult[["info"]][["directory"]],
                            qcResult[["info"]][["logFilename"]])
    plotF <- Sxslt::xsltApplyStyleSheet(plotFPath,
                                  xslStyles[["plotFunAndFileStyleSheet"]])
    logName <- logNameToHTML(plotFPath)
    saveXML(plotF$doc, file = logName)
    logName    
}

# --------------------------------------------------------------------
#
# report.qcCompareExprResult()
#
# --------------------------------------------------------------------
`report.qcCompareExprResult` <-
function(qcResult, xslStyles)
{
    compareExprPath <- file.path(qcResult[["info"]][["path"]],
                                 qcResult[["info"]][["logFilename"]])
    testPath <- file.path(qcResult[["testInfo"]][["directory"]],
                          qcResult[["testInfo"]][["logFilename"]])
    controlPath <- file.path(qcResult[["controlInfo"]][["directory"]],
                             qcResult[["controlInfo"]][["logFilename"]])

    compareExpr <- Sxslt::xsltApplyStyleSheet(compareExprPath,
                                       xslStyles[["compareExprStyleSheet"]])
    # No point recursing since we have the info we need and there's
    # only 2 files.
    testExpr <- Sxslt::xsltApplyStyleSheet(testPath,
                                    xslStyles[["plotExprStyleSheet"]])
    controlExpr <- Sxslt::xsltApplyStyleSheet(controlPath,
                                       xslStyles[["plotExprStyleSheet"]])
    logName <- logNameToHTML(compareExprPath)
    saveXML(compareExpr$doc, file = logName)
    saveXML(testExpr$doc, file = logNameToHTML(testPath))
    saveXML(controlExpr$doc, file = logNameToHTML(controlPath))
    logName
}

# --------------------------------------------------------------------
#
# report.qcCompareFileResult() and report.qcCompareFunResult()
#
# --------------------------------------------------------------------
`report.qcCompareFileResult` <- `report.qcCompareFunResult` <-
function(qcResult, xslStyles)
{
    lapply(qcResult[["results"]], report, xslStyles)

    compareFPath <- file.path(qcResult[["info"]][["path"]],
                              qcResult[["info"]][["logFilename"]])
    testPath <- qcResult[["info"]][["testLog"]]
    controlPath <- qcResult[["info"]][["controlLog"]]
    compareF <- Sxslt::xsltApplyStyleSheet(compareFPath,
                                xslStyles[["compareFunAndFileStyleSheet"]])
    # No point recursing, we have what we need.
    testF <- Sxslt::xsltApplyStyleSheet(testPath,
                                xslStyles[["plotFunAndFileStyleSheet"]])
    controlF <- Sxslt::xsltApplyStyleSheet(controlPath,
                                xslStyles[["plotFunAndFileStyleSheet"]])
    logName <- logNameToHTML(compareFPath)
    saveXML(compareF$doc, file = logName)
    saveXML(testF$doc, file = logNameToHTML(testPath))
    saveXML(controlF$doc, file = logNameToHTML(controlPath))
    logName
}

# --------------------------------------------------------------------
#
# report.default()
#
# --------------------------------------------------------------------
`report.default` <-
function(qcResult, xslStyles)
{
    stop(gettextf("comparing %s not yet implemented",
                  sQuote(class(qcResult))),
         domain=NA)
}

# Some miscellaneous functions used in this file and by the .xsl files
`logNameToHTML` <- function(logName) gsub("[.]xml$", ".html", logName)
`logToHTML` <- function(...) logNameToHTML(file.path(...))
`getCompareExprName` <- function(logWithPath) {
    strsplit(basename(logWithPath), "-compareExprLog.xml")[[1]]
}

