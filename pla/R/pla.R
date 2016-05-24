pla <- function(data,
                alpha            = 0.05,
                indexOfReference = 1,
                StdName          = sampleLabels[indexOfReference],
                sampleLabels     = "data",
                imputeMissing    = FALSE,
                dfAdjustment     = NA,
                dilutionRatio    = NA,
                factor           = NA,
                selectFun        = function (array) NULL,
                echoData         = TRUE,
                colors           = "default",
                projectTitle     = "",
                assayTitle       = "",
                design           = "", ...) {
    cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3",
                    "#0072B2", "#D55E00", "#CC79A7", "#F0E442")
    class <- "pla"
    if ((design == "") | (design == "fromData")) {
        if (any(names(data@labels) == "Design")) {
            design <- data@labels["Design"]
        }
        if (any(slotNames(data) == "design")) {
            design <- data@design
        }
    }
    design <- .string2design(design)
    if ((design == "lsd"))
        class <- "plaLSD"
    if ((design == "rbd"))
        class <- "plaRBD"
    if ((design == "crd") | (design == ""))
        class <- "plaCRD"
    if ((class(data) == "assayData") |
        (class(data) == "assayTable") |
        (class(data) == "assayFrame")) {
        frame <- data.frame()
        labels <- data@labels
        if (any(names(labels) == "DilutionRatio")) {
            dilutionRatio <- ifelse(is.na(dilutionRatio),
                                    ifelse(any(slotNames(data) ==
                                               "dilutionRatio"),
                                           data@dilutionRatio,
                                           as.numeric(labels["DilutionRatio"])),
                                    dilutionRatio)
        }
        dilutionRatio <- ifelse(is.na(dilutionRatio), 2.0, dilutionRatio)
        OK <- TRUE
        if (class(data) == "assayFrame") {
            OK <- .checkPlaFrame(data@tableRaw, design = design)
            ## if (design == "rbd")
            ##     warning(
            ##         "Random Block Designs expected to be entered as a table!")
            if (design == "lsd") {
                invalidLS <- .invalidLatinSquare(data@tableRaw)
                if (!is.null(invalidLS)) {
                    warning("Invalid Latin Square")
                    OK <- FALSE
                    print(invalidLS)
                    invisible(invalidLS)
                }
            }
            if (design == "rbd") {
                invalidRB <- .invalidBlockDesign(data@tableRaw)
                if (!is.null(invalidRB)){
                    warning("Invalid Block Design")
                    OK <- FALSE
                    print(invalidRB)
                    invisible(invalidRB)
                }
            }
            if (OK) {
                if (design == "lsd")
                    select <- c("Row", "Column", "Sample", "Dilution", "Response")
                else
                    select <- c("Replicate", "Sample", "Dilution", "Response")
                frame <- data@tableRaw
                frame <- data2assayFrame(frame[, select], dr = dilutionRatio,
                                         design = design)
                nMissing <- length(which(is.na(frame["Response"])))
            }
            if (design == "rbd")
                if (imputeMissing) {
                    warning(
                        "Imputation of Random Blocks not tested for data.frame!")
                   frame <- .imputeRBDframe(frame, nMissing = nMissing, ...)
                }
            if (design == "crd")
                if (imputeMissing)
                    frame <- .imputeCRD(frame, ...)
            if (design == "lsd")
               if (imputeMissing)
                   frame <- .imputeLatinSquare(frame, nMissing = nMissing, ...)
        } else {
            if (length(data@assay) > 0)
                table <- data@assay
            else
                table <- as.array(data, reduced = TRUE, selectFun = selectFun)
            nMissingTable <- length(which(is.na(table)))
            nMissing <- nMissingTable
            if (design == "lsd")
                    warning("Can you enter a Latin Square in this table?!?")
            if (design == "rbd")
                if (imputeMissing)
                    table <- .imputeRBD(table, nMissing = nMissingTable, ...)
            frame <- assayTable2frame(table, dr = dilutionRatio,
                                      echoData = echoData, ...)
            nMissingFrame <- length(which(is.na(frame["Response"])))
            if ((nMissingFrame != nMissingTable) & !(design == "rbd"))
                warning(
                    paste0("Number of missing values for table (", nMissingTable,
                           ") differs from number of missing vaules for frame (",
                           nMissingFrame, ")"))
            if (design == "crd")
                if (imputeMissing) {
                    frame <- .imputeCRD(frame, nMissing = nMissingTable, ...)
                }
            if (design == "lsd")
                if (imputeMissing) {
                    frame <- .imputeLatinSquare(frame,
                                                nMissing = nMissingTable, ...)
                }
        }
        if (OK)
            if (length(sampleLabels) == 1)
                if (sampleLabels == "data")
                    sampleLabels <- levels(unlist(frame["Sample"]))
        iOR <- match(StdName, sampleLabels)
        if ((indexOfReference == 1) & (indexOfReference != iOR))
            indexOfReference <- iOR
        if (colors == "default") {
            colors <- c("grey", cbbPalette[2:length(sampleLabels)])
        }
        if (any(names(labels) == "Factor") | any(names(labels) == "FactorA")) {
            if (length(factor) < 2)
                if (is.na(factor) & !is.nan(factor)) {
                    if (any(slotNames(data) == "factor"))
                        factor <- data@factor
                    else
                        factor <- as.numeric(labels["Factor"])
                    names(factor) <- ifelse(is.na(factor),
                                            "Rel. Est.", factor); ## ?!?
                }
        }
        if (length(factor) < 2)
            if (is.na(factor) & !is.nan(factor))
                factor <- 1.0
        if (any(names(labels) == "dfAdjustment") & is.na(dfAdjustment)) {
            dfAdjustment <- ifelse(is.na(dfAdjustment),
                                   ifelse(any(slotNames(data) == "dfAdjustment"),
                                          data@dfAdjustment,
                                          as.numeric(labels["dfAdjustment"])),
                                   dfAdjustment)
        }
        dfAdjustment <- ifelse(is.na(dfAdjustment), nMissing, dfAdjustment)
        new(class,
            data             = frame,
            alpha            = alpha,
            dfAdjustment     = dfAdjustment,
            dilutionRatio    = dilutionRatio,
            selectFun        = selectFun,
            factor           = factor,
            sampleLabels     = sampleLabels,
            indexOfReference = indexOfReference,
            StdName          = StdName,
            colors           = colors,
            projectTitle     = ifelse(projectTitle == "",
                data@projectTitle, projectTitle),
            assayTitle       = ifelse(assayTitle   == "",
                data@assayTitle,   assayTitle),
            imputeMissing    = imputeMissing,
            design           = design)
    }
    else if (class(data) == "data.frame") {
        dilutionRatio <- ifelse(is.na(dilutionRatio), 2.0, dilutionRatio)
        if (length(factor) < 2)
            if (is.na(factor) & !is.nan(factor))
                factor <- 1.0
        if (design == "lsd")
            select <- c("Row", "Column", "Sample", "Dilution", "Response")
        else
            select <- c("Replicate", "Sample", "Dilution", "Response")
        frame <- data.frame()
        OK <- .checkPlaFrame(data, design = design)
        ## if (design == "rbd")
        ##    warning("Random Block Designs expected to be entered as a table!")
        if (design == "lsd") {
            invalidLS <- .invalidLatinSquare(data[, select])
            if (!is.null(invalidLS)){
                warning("Invalid Latin Square")
                OK <- FALSE
                print(invalidLS)
                invisible(invalidLS)
            }
        }
        if (design == "rbd") {
            invalidRB <- .invalidBlockDesign(data[, select])
            if (!is.null(invalidRB)){
                warning("Invalid Block Design")
                OK <- FALSE
                print(invalidRB)
                invisible(invalidRB)
            }
        }
        if (OK)
            frame <- data2assayFrame(data[, select],
                                     dr = dilutionRatio, design = design)
        nMissing <- length(which(is.na(frame["Response"])))
        if (design == "rbd")
            if (imputeMissing) {
                warning(
                    "Imputation of Random Blocks not tested for data.frame!")
                frame <- .imputeRBDframe(frame, nMissing = nMissing, ...)
            }
        if (design == "crd")
            if (imputeMissing)
                frame <- .imputeCRD(frame, nMissing = nMissing, ...)
        if (design == "lsd")
            if (imputeMissing)
                frame <- .imputeLatinSquare(frame, nMissing = nMissing, ...)
        dfAdjustment <- ifelse(is.na(dfAdjustment), nMissing, dfAdjustment)
        if (OK)
            if (length(sampleLabels) == 1)
                if (sampleLabels == "data")
                    sampleLabels <- levels(unlist(frame["Sample"]))
        iOR <- match(StdName, sampleLabels)
        if ((indexOfReference == 1) & (indexOfReference != iOR))
            indexOfReference <- iOR
        if (colors == "default") {
            colors <- c("grey", cbbPalette[2:length(sampleLabels)])
        }
        new(class,
            data             = frame,
            alpha            = alpha,
            dfAdjustment     = dfAdjustment,
            dilutionRatio    = dilutionRatio,
            selectFun        = selectFun,
            factor           = factor,
            sampleLabels     = sampleLabels,
            indexOfReference = indexOfReference,
            StdName          = StdName,
            colors           = colors,
            projectTitle     = projectTitle,
            assayTitle       = assayTitle,
            imputeMissing    = imputeMissing,
            design           = design)
    }
            
}
