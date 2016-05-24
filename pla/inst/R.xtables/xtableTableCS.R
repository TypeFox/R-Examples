xtableTableCS <- function(imputedData,
                          data,
                          FUN = NULL,
                          label = "Data",
                          digits = 4,
                          response = "Response",
                          row = "Row",
                          replicate = "Replicate",
                          column = "Column",
                          sampleStep = "SampleStep",
                          sample = "Sample",
                          step = "Dilution",
                          ncol = 6,
                          doMeans = TRUE) {

    if (is.null(FUN)) {
        fun <- function(x) x
        if (!is.null(data))
            if (any(names(data@labels) == "INV")) {
                warning("Found function")
                fun <- data@inv
            }
    } else fun <- FUN

    if (!any(dimnames(imputedData)[[2]] == sampleStep)) {
        ## SampleStep <- paste0(as.character(unlist(imputedData[sample])), ":",
        ##                      as.character(unlist(imputedData[step])))
        SampleStep <- pasteSS(as.character(unlist(imputedData[sample])),
                              as.character(unlist(imputedData[step])))
        names(SampleStep) <- "SampleStep"
        imputedData <- cbind(imputedData, SampleStep = SampleStep)
    }

    strCaption <-
        paste0("\\textbf{Responses by Treatments (Sample \\& Dose) and Replicate}")
    ## D <- split(imputedData[, c(response)], imputedData[, sampleStep])
    D <- split(imputedData[, c(response)], imputedData[, replicate])
    E <- split(imputedData[, c(sampleStep)], imputedData[, replicate])
    table <- matrix(unlist(D), byrow = TRUE, nrow = length(D))
    table <- fun(table)
    nrows <- dim(table)[1]
    dimnames(table)[[1]] <- names(D)
    ncolumns <- dim(table)[2]
    dimnames(table)[[2]] <- E[[1]]
    ## table <- t(table)
    table <- rbind(table,
                   Mean =      apply(table, 2, mean, na.rm = TRUE),
                   SD   = sqrt(apply(table, 2, var , na.rm = TRUE)),
                   CV   = sqrt(apply(table, 2, var , na.rm = TRUE)) /
                   apply(table, 2, mean, na.rm = TRUE) * 100)
    table <- cbind(table,
                   Mean =      apply(table, 1, mean, na.rm = TRUE),
                   SD   = sqrt(apply(table, 1, var , na.rm = TRUE)),
                   CV   = sqrt(apply(table, 1, var , na.rm = TRUE)) /
                   apply(table, 1, mean, na.rm = TRUE) * 100)
    table[nrows + 1, ncolumns + 2] <- NA
    table[nrows + 1, ncolumns + 3] <- NA
    table[nrows + 2, ncolumns + 1] <- NA
    table[nrows + 2, ncolumns + 2] <- NA
    table[nrows + 2, ncolumns + 3] <- NA
    table[nrows + 3, ncolumns + 1] <- NA
    table[nrows + 3, ncolumns + 2] <- NA
    table[nrows + 3, ncolumns + 3] <- NA
    nms <- dimnames(table)[[2]]
    pos <- list(-1, dim(table)[1]-3, dim(table)[1])
    labelsTable <- labelsXfun(nms, pos)
    XT <- xtable(table, digits = c(117, rep(digits, dim(table)[2])),
                 caption = strCaption, label = "Data")
    align(XT) <- paste0(c("r|", rep("r", ncolumns), "|rrr|"), collapse = "")
    print(XT, size = "footnotesize", include.colnames = FALSE,
          hline.after = NULL, caption.placement = "top",
          sanitize.text.function = function(x) paste0(x),
          add.to.row = list(pos = pos, command = labelsTable))

    if (doMeans) {
        strCaption <- paste0("\\textbf{Means}")
        A <- split(imputedData[, c(response, step)], imputedData[, sample])
        B <- lapply(A, FUN = function(x) split(x[, response], x[, step]))
        C <- lapply(B, FUN = function(x) lapply(x, mean))
        table <- matrix(unlist(C), byrow = TRUE, nrow = length(A))
        nrows <- dim(table)[1]
        dimnames(table)[[1]] <- names(A)
        ncolumns <- dim(table)[2]
        dimnames(table)[[2]] <- names(B[[1]])
        nms <- dimnames(table)[[2]]
        pos <- list(-1, dim(table)[1]-3, dim(table)[1])
        labelsTable <- labelsXfun(nms, pos)
        XT <- xtable(table, digits = c(117, rep(2, dim(table)[2])),
                     caption = strCaption, label = "Data")
        print(XT, size = "footnotesize", include.colnames = FALSE,
              hline.after = NULL, caption.placement = "top",
              add.to.row = list(pos = pos, command = labelsTable))
    }
    
}
