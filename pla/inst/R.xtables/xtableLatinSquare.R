xtableLatinSquare <- function(data,
                              response = "Response",
                              row = "Row",
                              column = "Column",
                              sampleStep = "SampleStep",
                              sample = "Sample",
                              step = "Dilution",
                              doMeans = FALSE) {


    if (!any(dimnames(data)[[2]] == sampleStep)) {
        ## SampleStep <- paste0(as.character(unlist(data[sample])), ":",
        ##                     as.character(unlist(data[step])))
        SampleStep <- pasteSS(as.character(unlist(data[sample])),
                              as.character(unlist(data[step])))
        names(SampleStep) <- "SampleStep"
        data <- cbind(data, SampleStep = SampleStep)
    }

    strCaption <-
        paste0("\\textbf{Treatments (Sample \\& Dose) by Row and Column}")
    As <- split(data[, c(sample, step, column)], data[, row])
    ncol <- length(As)
    Bs <- unlist(lapply(As, FUN = function(x)
                        split(pasteSS(as.character(x[, sample]),
                                      as.character(x[, step])), x[, column])))
    tableLabels <- matrix(Bs, ncol = ncol)
    nrows <- dim(tableLabels)[1]
    dimnames(tableLabels)[[1]] <- paste("", 1:dim(tableLabels)[1], sep = " ")
    ncolumns <- dim(tableLabels)[2]
    dimnames(tableLabels)[[2]] <- paste("", 1:dim(tableLabels)[2], sep = " ")
    nms <- dimnames(tableLabels)[[2]]
    pos <- list(-1, dim(tableLabels)[1])
    labelsTable <- labelsXfun(nms, pos)
    XT <- xtable(tableLabels, caption = strCaption, label = "Data")
    print(XT, size = "footnotesize", include.colnames = FALSE,
          hline.after = NULL, caption.placement = "top",
          sanitize.text.function = function(x) paste0(x),
          add.to.row = list(pos = pos, command = labelsTable))

    A <- split(data[, c(response, column)], data[, row])
    B <- unlist(lapply(A, FUN = function(x) split(x[, response], x[, column])))
    table <- matrix(B, byrow = TRUE, ncol = ncol)
    nrows <- dim(table)[1]
    dimnames(table)[[1]] <- paste("", 1:dim(table)[1], sep = " ")
    ncolumns <- dim(table)[2]
    dimnames(table)[[2]] <- paste("", 1:dim(table)[2], sep = " ")
    table <- rbind(table,
                   Mean =    apply(table, 2, mean, na.rm = TRUE),
                   SD = sqrt(apply(table, 2, var,  na.rm = TRUE)),
                   CV = sqrt(apply(table, 2, var,  na.rm = TRUE)) /
                   apply(table, 2, mean, na.rm = TRUE) * 100)

    table <- cbind(table,
                   Mean =    apply(table, 1, mean, na.rm = TRUE),
                   SD = sqrt(apply(table, 1, var,  na.rm = TRUE)),
                   CV = sqrt(apply(table, 1, var,  na.rm = TRUE)) /
                   apply(table, 1, mean, na.rm = TRUE) * 100)
    table[nrows + 1, ncolumns + 2] <- NA
    table[nrows + 1, ncolumns + 3] <- NA
    table[nrows + 2, ncolumns + 1] <- NA
    table[nrows + 2, ncolumns + 2] <- NA
    table[nrows + 2, ncolumns + 3] <- NA
    table[nrows + 3, ncolumns + 1] <- NA
    table[nrows + 3, ncolumns + 2] <- NA
    table[nrows + 3, ncolumns + 3] <- NA

    strCaption <- paste0("\\textbf{Responses Row and Column (Replication)}")
    nms <- dimnames(table)[[2]]
    pos <- list(-1, dim(table)[1]-3, dim(table)[1])
    labelsTable <- labelsXfun(nms, pos)
    XT <- xtable(table, digits = c(117, rep(2, dim(table)[2])),
                 caption = strCaption, label = "Data")
    align(XT) <- paste0(c("r|", rep("r", nrows), "|rrr|"), collapse = "")
    print(XT,size = "footnotesize", include.colnames = FALSE,
          hline.after = NULL, caption.placement = "top",
          add.to.row = list(pos = pos, command = labelsTable))

    strCaption <-
        paste0("\\textbf{Responses by Treatments (Sample \\& Dose) and Repetition}")
    D <- split(data[, c(response)], data[, sampleStep])
    table <- matrix(unlist(D), byrow = TRUE, nrow = length(D))
    nrows <- dim(table)[1]
    dimnames(table)[[1]] <- names(D)
    ncolumns <- dim(table)[2]
    dimnames(table)[[2]] <- 1:ncolumns
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
    XT <- xtable(table, digits = c(117, rep(2, dim(table)[2])),
                 caption = strCaption, label = "Data")
    align(XT) <- paste0(c("r|", rep("r", nrows), "|rrr|"), collapse = "")
    print(XT,size = "footnotesize", include.colnames = FALSE,
          hline.after = NULL, caption.placement = "top",
          sanitize.text.function = function(x) paste0(x),
          add.to.row = list(pos = pos, command = labelsTable))

    if (doMeans) {
        strCaption <- paste0("\\textbf{Means}")
        A <- split(data[, c(response, step)], data[, sample])
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
