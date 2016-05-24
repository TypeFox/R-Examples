
## Method: 'as.table':
## Method: 'as.array':
## Method: 'as.data.frame':

## Method: 'show':      Display the object, by printing, plotting or
##                      whatever suits its class
## Method: 'print':     Prints its argument and returns it invisibly
##                      (via invisible(x))

setMethod("as.table",
    signature(x = "assayFrame"),
    function (x, which = 1:length(factorNames), type = "counts",
              response = responseName,
              reduced = FALSE, selectFun = function (array) NULL, ...)
    {
        if (length(x@internals) > 0) {
            indexHeader <- x@internals$description$indexHeader
            responseName <- names(indexHeader)[1]
            exclude <- names(indexHeader) == response |
                names(indexHeader) == "Label"
            factorNames <- unlist(indexHeader[!exclude])
            ## if (type == "counts") {
            ##     table(x@tableRaw[, factorNames[rev(which)]])
            ## } else {
            ## }
            as.array(x, which = which, type = type,
                     reduced = reduced, selectFun = selectFun, ...)
        } else
            as.array(x, which = which, type = type,
                     reduced = reduced, selectFun = selectFun, ...)
    }
)

setMethod("as.table",
    signature(x = "assayTable"),
    function (x, ...)
    {
        ifelse(is.na(as.array(x, ...)), 0, 1)
    }
)

##        .f <- function(X, dmsT, response = "Response") {
##            factorname <- names(dmsT)[1]
##            Result <- NULL
##            for (level in dmsT[[1]]) {
##                select <- X[, factorname] == level
##                if (length(dmsT) > 1)
##                    Result <- append(Result, f(X[select,],
##                                               dmsT[-1], response))
##                else
##                    Result <- append(Result, FUN(X[select, response]))
##            }
##            return(Result)
##        }

setMethod("as.array",
    signature(x = "data.frame"),
    function (x, whichFactors = 1:length(factorNames),
              whichResponse = length(responseNames),
              type = "values", FUN = function(x) mean(x),
              response = responseNames[whichResponse],
              responseNames = names(which(lapply(x, class) == "numeric")),
              factorNames = names(which(lapply(x, class) == "factor")),
              reduced = FALSE, selectFun = function (array) NULL,
              ...)
    {
        f <- function(X, dmsT, response = "Response",
                      factorname = names(dmsT)[1])
            return(unlist(lapply(dmsT[[1]],
                                 FUN = function(level,
                                     select = X[, factorname] == level)
                                     return(ifelse(length(dmsT) > 1,
                                                   list(f(X[select,],
                                                          dmsT[-1], response)),
                                                   FUN(X[select, response]))))))
        if ((class(whichFactors) == "integer") |
            (class(whichFactors) == "numeric"))
            selectNames <- factorNames[whichFactors]
        else
            selectNames <- whichFactors
        ## Hack for 'x' is from 'assayTable2frame' :
        if (length(factorNames) == length(which(lapply(x, class) == "factor")))
            if (all(factorNames == names(which(lapply(x, class) == "factor"))))
                if ((length(selectNames) == 2) &
                    any("Dilution" == dimnames(x)[[2]]))
                    if (all(selectNames == c("SampleStep", "Sample")))
                        selectNames <- c("Sample", "Dilution")
        Array <- table(x[, selectNames])
        dimNames <- dimnames(Array)
        if (type != "counts") {
            Responces <- f(x, rev(dimNames), response)
            Array <- array(Responces, dim(Array), dimnames = dimNames)
        }
        if (!is.null(body(selectFun)))
            Array <- selectFun(Array)
        if (reduced)
            Array <- .reduceDimension(
                Array, dimNames, selectFun = selectFun)
        return(Array)
    }
)

setMethod("as.array",
    signature(x = "assayFrame"),
    function (x, which = 1:length(factorNames),
              type = "values", FUN = function(x) mean(x),
              response = responseName,
              reduced = FALSE, selectFun = function (array) NULL,
              ...)
    {
        if (length(x@assay) > 0)
            return(x@assay)
        else {
            if (length(x@internals) > 0) {
                indexHeader <- x@internals$description$indexHeader
                responseName <- names(indexHeader)[1]
                exclude <- names(indexHeader) == response |
                    names(indexHeader) == "Label"
                factorNames <- unlist(indexHeader[!exclude])
                Array <- as.array(x@tableRaw, whichFactors = which,
                                  whichResponse = 1, type = type,
                                  FUN = FUN,
                                  response = responseName,
                                  responseNames = responseName,
                                  factorNames = factorNames,
                                  reduced = reduced,
                                  selectFun = selectFun, ...)
                ## ## Not relevant here, as the 'responses' are as the are:
                ## fun <- x@fun
                ## log <- x@log
                ## Array <- fun(Array)
                ## if (log == TRUE)
                ##     Array <- log(Array)
                ## else if (log != FALSE)
                ##     Array <- log(Array) / log(log)
                ## ## Alternativ: Then 'fun' and 'log' should also be applied on
                ## ##             'dataframe["Dilution"]' in 'data2assayFrame'!
                return(Array)
            } else
                return(as.array(x@tableRaw, type = type, FUN = FUN,
                                ## which = which, response = response,
                                reduced = reduced, selectFun = selectFun, ...))
        }
    }
)

setMethod("as.array",
    signature(x = "assayTable"),
    function (x, reduced = FALSE, selectFun = function (array) NULL, ...)
    {
        if (length(x@assay) > 0)
            return(x@assay)
        else {
            internals <- x@internals
            AD <- internals$description
            TableChar <- x@tableRaw
            fun <- x@fun
            log <- x@log
            if (!is.null(dim(TableChar)))
                if (dim(TableChar)[1] > 0) {
                    Array <- .apermAssay(
                        TableChar, fun, log, AD$inversPerm,
                        AD$Dimension, AD$orderedNames,
                        combinedTreatment = AD$combinedTreatment)
                    if (!is.null(body(selectFun)))
                        Array <- selectFun(Array)
                    if (reduced)
                        Array <- .reduceDimension(
                            Array, AD$orderedNames[AD$inversPerm],
                            selectFun = selectFun)
                    return(Array)
                }
        }
    }
)

setMethod("as.data.frame",
    signature(x = "assayFrame"),
    function (x, ...)
    {
        return(x@tableRaw)
    }
)

setMethod("as.data.frame",
    signature(x = "assayTable"),
    function (x, row.names = NULL, optional = FALSE, ...)
    {
        return(assayTable2frame(as.array(x, reduced = TRUE, ...),
                                ## dr = x@dilutionRatio,
                                ...))
    }
)

setMethod("show",
    signature(object = "assayFrame"),
    function (object)
    {
        ih <- object@internals$description$indexHeader
        show(str(ih[lapply(ih, length) != 0]))
        printOption <- function(label, txt)
            if (length(txt) > 0)
                if (any(txt != ""))
                    cat(paste0(label, txt, "\n"))
        printOption("Project:       ", object@projectTitle  )
        printOption("Assay:         ", object@assayTitle    )
        printOption("Description:   ", object@description   )
        printOption("Comment:       ", object@comment       )
        printOption("Resume:        ", object@resume        )
        printOption("Date:          ", object@date          )
        printOption("Operator:      ", object@operator      )
        printOption("Model:         ", object@model         )
        printOption("Design:        ", object@design        )
        printOption("Factor(s):     ",
                    paste(prettyNum(object@factor), collapse = ", "))
        printOption("Dilution ratio:", object@dilutionRatio )
        printOption("Adjustment, DF:", object@dfAdjustment  )
        Array <- as.array(object)
        show(ftable(Array, col.vars = length(dim(Array)) + c(-1, 0)))

    }
)

setMethod("show",
    signature(object = "assayTable"),
    function (object)
    {
        show(str(object@internals$description$orderedNames))
        printOption <- function(label, txt)
            if (length(txt) > 0)
                if (any(txt != ""))
                    cat(paste0(label, txt, "\n"))
        printOption("Project:       ", object@projectTitle  )
        printOption("Assay:         ", object@assayTitle    )
        printOption("Description:   ", object@description   )
        printOption("Comment:       ", object@comment       )
        printOption("Resume:        ", object@resume        )
        printOption("Date:          ", object@date          )
        printOption("Operator:      ", object@operator      )
        printOption("Model:         ", object@model         )
        printOption("Design:        ", object@design        )
        printOption("Factor(s):     ",
                    paste(prettyNum(object@factor), collapse = ", "))
        printOption("Dilution ratio:", object@dilutionRatio )
        printOption("Adjustment, DF:", object@dfAdjustment  )
        ## Array <- as.array(object)
        ## show(ftable(Array, col.vars = length(dim(Array)) + c(-1, 0)))
        Array <- as.array(object, reduced = TRUE)
        show(ftable(Array))
    }
)

setMethod("print",
    signature(x = "assayFrame"),
    function (x, ...)
    {
        print(x@tableRaw)
    }
)

setMethod("print",
    signature(x = "assayTable"),
    function (x, oldForm = FALSE, ...)
    {

        if (!oldForm)
            if (length(x@assay) > 0)
                print(x@assay)
            else
                print(x@tableRaw)
        else {
            internals <- x@internals
            TableChar <- x@tableRaw
            fun <- x@fun
            log <- x@log
            if (!is.null(dim(TableChar)))
                if (dim(TableChar)[1] > 0) {
                    Table <- .transposeTable(
                        TableChar, internals$firstName,
                        internals$firstColumn,
                        internals$columnNames, fun, log)
                    print(Table)
                }
        }
    }
)
