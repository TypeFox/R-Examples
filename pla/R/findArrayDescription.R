.findArrayDescription <-
    function (Data,
              Plates,
              TableChar,
              HeaderMatch,
              FactorsMatch,
              indexMatch,
              rows = "",
              columns = "",
              replicateNames = c("Repetition", "Replicate", "Block",
                  "Plate", "Animal", "Group", "Litter"),
              doseNames = c("Dilution", "Dose"),
              sampleNames = c("Smp", "Sample"),
              otherNames = c("Outer", "Inner", "Assay"),
              combinedTreatment = FALSE)
{
    Error <- 0
    if (length(HeaderMatch[FactorsMatch != "Label"]) > 0) {
        strSplitHeaders   <- strsplit(dimnames(TableChar)[[2]], "[.]")
        lengthItems <- lapply(strSplitHeaders, length)
        if (diff(range(unlist(lengthItems))) != 0) {
            ## Correct head 'read.table' adds "_1" to
            ## first replicate of column label:
            altSplitHeaders <-
                lapply(strSplitHeaders, FUN = function(x)
                       if (length(x) == 1) x else
                       if (!is.na(match(x[1], HeaderMatch)) &
                           (x[2] == "1")) x[1] else x)
            lengthAltItems <- lapply(altSplitHeaders, length)
            if (diff(range(unlist(lengthAltItems))) == 0) {
                strSplitHeaders <- altSplitHeaders
                lengthItems <- lengthAltItems
            }
        }
        if (all(unlist(lengthItems) == 1)) {
            if (all(unlist(lapply(strSplitHeaders, FUN = nchar)) == 1))
                matrixHeaders <- matrix(
                    ## If only one character, then an 'empty' columns:
                    unlist(lapply(strSplitHeaders, FUN = function(x)
                        c(substr(x, 1, 1)))),
                    byrow = TRUE, nrow = length(dimnames(TableChar)[[2]]))
            else
                matrixHeaders <- matrix(
                    unlist(lapply(strSplitHeaders, FUN = function(x)
                        c(substr(x, 1, 1), substr(x, 2, nchar(x))))),
                    byrow = TRUE, nrow = length(dimnames(TableChar)[[2]]))
        }
        else {
            if (diff(range(unlist(lengthItems))) != 0) {
                warning("Variating number of items in headers",
                        call. = FALSE)
                strSplitHeaders <- lapply(strSplitHeaders,
                                          FUN = function(x)
                                              c(x, rep(" ",
                                                       max(unlist(lengthItems))
                                                       - length(x))))
                warning(paste("Please use identical patterns",
                              "of special characters \n",
                              " in factor level specifications in the header",
                              collapse = ""),
                        call. = FALSE)
                warning(paste(paste(dimnames(TableChar)[[2]], collapse = ", "),
                              "replaced by \n ",
                              paste(lapply(strSplitHeaders,
                                           FUN = function(x)
                                               paste(x, collapse = ".")),
                                    collapse = ", "),
                              collapse = " "),
                        call. = FALSE)
            }
            matrixHeaders <- matrix(unlist(strSplitHeaders), byrow = TRUE,
                                    nrow = length(strSplitHeaders))
        }
        nHeaders          <- dim(matrixHeaders)[[1]]
        namesHeaders      <- apply(matrixHeaders, 2, FUN = function(x)
            paste0(unique(x), collapse = "/"))
        describeHeader    <- apply(matrixHeaders, 2, FUN = function(x)
            .countFactorLevels(x, warnings = FALSE)[[1]])
        if (dim(matrixHeaders)[2] > 2)
            if (!all(describeHeader["prod",] == dim(matrixHeaders)[1]) | TRUE) {
                for (i in (dim(matrixHeaders)[2]):3) {
                    t1 <- .is.digits(matrixHeaders[, i  ])
                    t2 <- .is.digits(matrixHeaders[, i-1])
                    v  <- unique(matrixHeaders[, i-2])
                    t3 <- FALSE
                    if (length(v) == 1)
                        if (is.character(v))
                            if ((tolower(v) == "d") | (tolower(v) == "dose"))
                                t3 <- TRUE
                    if (t1 & t2 & t3) {
                        dose <- paste(as.character(matrixHeaders[, i-1]),
                                      as.character(matrixHeaders[, i]),
                                      sep = ".")
                        matrixHeaders[, i-1] <- dose
                        matrixHeaders <- matrixHeaders[, -i]
                    }
                }
                namesHeaders      <- apply(matrixHeaders, 2, FUN = function(x)
                    paste0(unique(x), collapse = "/"))
                describeHeader    <- apply(matrixHeaders, 2, FUN = function(x)
                    .countFactorLevels(x)[[1]])
            }
        dimnames(describeHeader)[[2]] <- namesHeaders
        matrixRows        <- Data[Plates, indexMatch[FactorsMatch != "Label"]]
        if (!is.null(dim(matrixRows)))
            if (dim(matrixRows)[1] < 1)
                warning("No data (levels of first factor invalid)")
        if (is.null(dim(matrixRows))) {
            matrixRows    <- unlist(strsplit(unlist(strsplit(
                as.character(matrixRows), "=")), "-"))
            matrixRows    <- matrix(matrixRows, byrow = TRUE,
                                    nrow = length(which(Plates)))
        }
        else
            matrixRows    <- t(apply(matrixRows, 1, FUN = function(x)
                unlist(strsplit(unlist(strsplit(x, "=")), "-"))))
        nRows             <- dim(matrixRows)[[1]]
        namesRows         <- apply(matrixRows, 2, FUN = function(x)
            paste0(unique(x), collapse = "/"))
        describeRows      <- apply(matrixRows, 2, FUN = function(x)
            .countFactorLevels(x)[[1]])
        dimnames(describeRows)[[2]] <- namesRows
        describeRows   <- .patchDescription(describeRows,   nRows,    "COLUMNS")
        describeHeader <- .patchDescription(describeHeader, nHeaders, "HEADER")
        bolRow <- describeRows["inter", ] == 1
        idxRow <- (1:dim(describeRows)[[2]])[bolRow]
        bolHed <- describeHeader["inter", ] == 1
        idxHed <- (1:dim(describeHeader)[[2]])[bolHed]
        Position <- c(idxRow+1, - idxHed - 1)
        names(Position) <- c(dimnames(describeRows)[[2]][bolRow],
                             dimnames(describeHeader)[[2]][bolHed])
        if (length(idxHed) > 0) {
            matrixHeadersRest <- matrixHeaders[, -c(idxHed, idxHed+1)]
            idxHedRest <- (1:dim(matrixHeaders)[[2]])[-c(idxHed, idxHed+1)]
        } else {
            matrixHeadersRest <- matrixHeaders
            idxHedRest <- (1:dim(matrixHeaders)[[2]])
        }
        if (length(idxRow) > 0) {
            matrixRowsRest <- matrixRows[, -c(idxRow, idxRow+1)]
            idxRowRest <- (1:dim(matrixRows)[[2]])[-c(idxRow, idxRow+1)]
        } else {
            matrixRowsRest <- matrixRows
            idxRowRest <- (1:dim(matrixRows)[[2]])
        }
        mcm <- function(x, table = names(Position)) {
            m <- charmatch(substring(tolower(x), 1, 3),
                           substring(tolower(table), 1, 3))
            m <- c(m, charmatch(substring(tolower(x), 1, 1),
                                substring(tolower(table), 1, 1)))
            m <- m[!is.na(m)]
            return(unique(sort(m)))
        }
        indexPosition <- list(
            Inner     = mcm(c("Inner")),
            Sample    = mcm(sampleNames),
            Dose      = mcm(doseNames),
            Replicate = mcm(replicateNames),
            Outer     = mcm(c("Outer")),
            Assay     = mcm(c("Assay")),
            Label     = mcm(c("Label")))
        ## Patch for extra factors is just theses four plus the below two lines:
        restNames <- otherNames[is.na(match(otherNames,
                                            c("Outer", "Inner", "Assay")))]
        if (length(restNames) > 0)
            for (i in 1:length(restNames)) {
                indexPosition <- append(indexPosition, mcm(restNames[i]))
                names(indexPosition)[length(indexPosition)] <- restNames[i]
            }
        missingFactors <- lapply(indexPosition[2:4], length) == 0
        if (length(which(missingFactors)) > 0) {
            if (length(which(missingFactors)) == 1) {
                if ((length(idxHedRest) + length(idxRowRest)) == 1) {
                    indexPosition[which(missingFactors)] <- length(Position) + 1
                    if (length(idxHedRest) == 1)
                        Position <- c(Position, -idxHedRest)
                    else
                        Position <- c(Position, idxRowRest)
                    names(Position) <- c(dimnames(describeRows)[[2]][bolRow],
                                         dimnames(describeHeader)[[2]][bolHed],
                                         names(which(missingFactors)))

                }  else {
                    Error <- Error + 1
                    warning(paste("More than one choice or none choices ",
                                  "for factor", names(which(missingFactors))),
                            call. = FALSE)
                }
            }
            if (length(which(missingFactors)) == 2) {
                cat(paste("Note: Two factors without labels: ",
                          paste(names(which(missingFactors)),
                                collapse = ", "), ".\n"))
                if ((length(idxHedRest) + length(idxRowRest)) == 2) {
                    if (length(idxHedRest) == 2)
                        {
                            indexPosition[which(missingFactors)[1]] <-
                                length(Position) + 1
                            indexPosition[which(missingFactors)[2]] <-
                                length(Position) + 2
                            if (describeHeader["class",idxHedRest][1] >
                                describeHeader["class",idxHedRest][2])
                                Position <- c(Position, -idxHedRest)
                            else
                                Position <- c(Position, -rev(idxHedRest))
                            names(Position) <- c(
                                dimnames(describeRows)[[2]][bolRow],
                                dimnames(describeHeader)[[2]][bolHed],
                                names(which(missingFactors)))

                        } else
                            warning(paste("Please use labels for rows: ",
                                          paste(names(which(missingFactors)),
                                                collapse = ", ")),
                                    call. = FALSE)
                } else {
                    Error <- Error + 2
                    warning(paste("More choices than missing factors:",
                                  paste(names(which(missingFactors)),
                                        collapse = ", ")), call. = FALSE)
                }
            }
            if (length(which(missingFactors)) == 3) {
                Error <- Error + 4
                warning(paste("To many factors (All!) without labels: ",
                              paste(names(which(missingFactors)),
                                    collapse = ", ")), call. = FALSE)
            }
        }
        posRow <- Position[Position > 0]
        valRow <- describeRows["before", posRow] -
            1 / describeRows["after", posRow]
        dimRow <- describeRows["inter", posRow]
        nmsRow <- matrixRows[, posRow]
        posHed <- Position[Position < 0]
        valHed <- describeHeader["before", -posHed] -
            1 / describeHeader["after", -posHed]
        dimHed <- describeHeader["inter", -posHed]
        nmsHed <- matrixHeaders[, -posHed]
        orderedPosition <- c(posRow[order(-valRow)], posHed[order(-valHed)])
        ## orderedNames <- list(nmsRow[, order(-valRow)],
        ##                      nmsHed[, order(-valHed)])
        if (prod(dimRow) != nRows)
            warning("Some factors have problems with nesting in columns",
                    call. = FALSE)
        if (prod(dimHed) != nHeaders)
            warning("Some factors have problems with nesting in header",
                    call. = FALSE)
        toList <- function(X) lapply(apply(X, 2,
                                           FUN = function(x) list(x)),
                                     FUN = function(x) x[[1]])
        if (length(posRow) > 1)
            orderedNamesR     <- toList(nmsRow[, order(-valRow)])
        else
            orderedNamesR     <- list(nmsRow)
        if (length(posHed) > 1)
            orderedNamesC     <- toList(nmsHed[, order(-valHed)])
        else
            orderedNamesC     <- list(nmsHed)
        orderedNames        <- append(orderedNamesR, orderedNamesC)

        Dimensions0         <- c(prod(dimRow[order(-valRow)]),
                                 prod(dimHed[order(-valHed)]),
                                 dimRow[order(-valRow)], dimHed[order(-valHed)])
        Dimensions          <- c(nRows, nHeaders,
                                 dimRow[order(-valRow)], dimHed[order(-valHed)])
        names(Dimensions0)  <- c("ROW", "COLUMN", names(orderedPosition))
        names(Dimensions)   <- c("ROW", "COLUMN", names(orderedPosition))
        names(orderedNames) <- names(orderedPosition)
        nmsOrdPos <- substr(names(orderedPosition), 1, 1)
        ## 'Inner' then should be eliminated by 'fun' in 'apermAssay':
        nmsOrdPos <- nmsOrdPos[nmsOrdPos != "I"]
        ## We want 'Replicates', 'Doses', and 'Samples' (plus 'Inner')
        ## to be 'last':
        perm <- ifelse(nmsOrdPos == "O", 0,
                       ifelse(nmsOrdPos == "A", 1,
                              ifelse(nmsOrdPos == "R", 2,
                                     ifelse(nmsOrdPos == "D", 3,
                                            ifelse(nmsOrdPos == "S", 4,
                                                   ifelse(nmsOrdPos == "I", 5,
                                                          NA))))))
        userFactors <- is.na(perm)
        ## Patch for extra factors is just theses two lines
        ## after the above four lines:
        perm[userFactors] <- 1 + (1:length(which(userFactors))) /
            (1 + length(which(userFactors)))
        names(perm) <- nmsOrdPos
        inversPerm <- order(perm)
        cmpNames <- function(namesArgs, namesData) {
            namesArgs[namesArgs == "Dilutions"] <- "Dose"
            for (name in replicateNames)
                namesArgs[namesArgs == name] <- "Replicate"
            for (name in doseNames)
                namesArgs[namesArgs == name] <- "Dose"
            for (name in sampleNames)
                namesArgs[namesArgs == name] <- "Sample"
            for (name in replicateNames)
                namesData[namesData == name] <- "Replicate"
            for (name in doseNames)
                namesData[namesData == name] <- "Dose"
            for (name in sampleNames)
                namesData[namesData == name] <- "Sample"
            if (length(namesArgs) != length(namesData))
                warning(paste(namesArgs, " / ",
                              paste(namesData, collapse = "/")),
                        call. = FALSE)
            else
                for (i in 1:length(namesArgs))
                    if (is.na(charmatch(namesArgs[i], namesData[i])) &
                        is.na(charmatch(namesData[i], namesArgs[i])))
                        warning(paste(namesArgs[i], " / ", namesData[i]),
                                call. = FALSE) }

        if (rows != "") {
            namesArgs <- unlist(strsplit(unlist(strsplit(rows, "/")), "&"))
            cmpNames(gsub(" ", "", namesArgs), names(posRow[order(-valRow)]))
        }
        if (columns != "") {
            namesArgs <- unlist(strsplit(unlist(strsplit(columns, "/")), "&"))
            cmpNames(gsub(" ", "", namesArgs), names(posHed[order(-valHed)]))
        }
        result <- list(matrixHeaders = matrixHeaders,
                       matrixRows = matrixRows,
                       combinedTreatment = combinedTreatment,
                       orderedPosition = orderedPosition,
                       orderedNames = orderedNames,
                       Dimensions = Dimensions,
                       inversPerm = inversPerm,
                       describeHeader = describeHeader,
                       describeRows = describeRows,
                       Error = Error)
    }
    else
        result <- NULL
    return(result)
}

.patchDescription <- function(Description, product, label = "COLUMNS") {
    if (prod(Description["inter", ]) != product) {
        ## warning(paste("Some factors have problems with nesting in", label),
        ##        call. = FALSE)
        cat(paste("Some factors have problems with nesting in", label, "\n"))
        labels <- Description["before", ] == 1 &
            Description["inter", ]  == 1 &
                Description["after", ]  == product
        invalid <- Description["before", ] == 1 &
            Description["inter", ]  == product &
                Description["after", ]  == 1
        rest  <- !(labels | invalid)
        if ((length(which(invalid)) == 1) & (length(which(rest)) == 1)) {
            Description["before", invalid] <- Description["inter",  rest]
            Description["inter",  invalid] <- Description["after",  rest]
            ## warning(paste("Trying to fix", names(which(invalid)), "by",
            ##              names(which(rest)), collapse = ""),
            ##        call. = FALSE)
            cat(paste("Trying to fix", names(which(invalid)), "by",
                      names(which(rest)), "\n", collapse = ""))
        }
    }
    return(Description)
}

## .my.as.numeric <- function(x) try(x <- as.numeric(x), silent = TRUE)
