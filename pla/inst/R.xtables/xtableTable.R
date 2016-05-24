xtableTable <- function(data,
                        digits = 4,
                        tabularx = FALSE,
                        verb = FALSE,
                        label = "Data",
                        INV = NULL,
                        sampleLabels = NULL,
                        transpose = FALSE,
                        selectFun = function (array) NULL) {

    if (is.null(INV)) {
        fun <- function(x) x
        if (any(names(data@labels) == "INV")) {
            fun <- data@inv
        }
    } else fun <- INV
    if (length(data@assay) > 0)
        TABLE <- data@assay
    else
        TABLE <- as.array(data, reduced = TRUE, selectFun = selectFun)
    table <- fun(TABLE)
    Labels1 <- updateLabels(table, 1)
    if (!is.null(Labels1))
        dimnames(table)[[1]] <- Labels1
    Labels2 <- updateLabels(table, 2)
    if (!is.null(Labels2))
        dimnames(table)[[2]] <- Labels2
    nrows <- dim(table)[1]
    ncols <- dim(table)[2]
    nSamples <- length(sampleLabels)
    nSteps <- ncols / nSamples

    extraSum <- FALSE
    if (data@design == "rbd" | data@design == "Randomised block")
        if (any(is.na(TABLE)))
            extraSum <- TRUE
    Names <- dimnames(TABLE)[[2]]
    tNames <- substr(Names, 1, 1)
    sNames <- unique(tNames)
    first <- match(sNames, tNames)
    if (max(first) > nSamples) {
        nLines <- nSamples-1
        nRows <- nSteps
    } else {
        nLines <- nSteps-1
        nRows <- nSamples
    }
    nr <- nrow(table)
    if (is.null(dimnames(table)[[1]]))
        dimnames(table)[[1]] <- 1:dim(table)[1]
    if (extraSum)
    table <- rbind(table,
                   Sum =     apply(table, 2,  sum, na.rm = TRUE))
    table <- rbind(table,
                   Mean =    apply(table[1:nr, ], 2, mean, na.rm = TRUE),
                   SD = sqrt(apply(table[1:nr, ], 2,  var, na.rm = TRUE)),
                   CV = sqrt(apply(table[1:nr, ], 2,  var, na.rm = TRUE)) /
                             apply(table[1:nr, ], 2, mean, na.rm = TRUE) * 100)
    if (extraSum) {
        table <- cbind(table,
                       Sum =     apply(table, 1,  sum, na.rm = TRUE),
                       Mean =    apply(table, 1, mean, na.rm = TRUE),
                       SD = sqrt(apply(table, 1, var,  na.rm = TRUE)),
                       CV = sqrt(apply(table, 1, var,  na.rm = TRUE)) /
                                 apply(table, 1, mean, na.rm = TRUE) * 100)
        table[nrows + 1, ncols + 2] <- NA
        table[nrows + 1, ncols + 3] <- NA
        table[nrows + 1, ncols + 4] <- NA
        table[nrows + 2, ncols + 1] <- NA
        table[nrows + 2, ncols + 2] <- NA
        table[nrows + 2, ncols + 3] <- NA
        table[nrows + 2, ncols + 4] <- NA
        table[nrows + 3, ncols + 1] <- NA
        table[nrows + 3, ncols + 2] <- NA
        table[nrows + 3, ncols + 3] <- NA
        table[nrows + 3, ncols + 4] <- NA
        table[nrows + 4, ncols + 1] <- NA
        table[nrows + 4, ncols + 2] <- NA
        table[nrows + 4, ncols + 3] <- NA
        table[nrows + 4, ncols + 4] <- NA
    }

    xtableHead(data, tabularx = tabularx, verb = verb) 
    
    dotable <- function(table, transpose = FALSE) {
        strCaption <- paste0("\\textbf{", label, "}")
        if (digits < 1) {
            table <- round(table, digits = - digits)
            ## http://jason.bryer.org/posts/2013-01-04/xtable_with_aligned_decimals.html
            XT <- xtable(matrix(prettyNum(table), ncol = ncol(table)),
                         caption = strCaption, label = label)
        } else {
            XT <- xtable(table, digits = c(117, rep(digits, dim(table)[2])),
                         caption = strCaption, label = label)
        }
        ## Does not work with the above 'updateLabels':
        ## nms <- toTex(dimnames(table)[[2]])
        nms <- (dimnames(table)[[2]])
        if (transpose) {
            if (digits == 0)
                align(XT) <- paste0(c("r|", rep("r", nrows), "|rrr|"), collapse = "")
            else
                align(XT) <- paste0(c("r|", rep("l", nrows), "|rrr|"), collapse = "")
            pos <- as.list(c(-1, nRows * (1:nLines), dim(table)[1]))
            labelsTable <- labelsXfun(nms, pos)
        } else {
            ps <- 3
            if (extraSum) {
                align(XT) <- paste0(c("r|", rep("r", ncols), "|rrrr|"), collapse = "")
                ps <- 4
            }
            pos <- list(-1, dim(table)[1]-ps, dim(table)[1])
            labelsTable <- labelsXfun(nms, pos)
        }
        print(XT, size = "footnotesize", include.colnames = FALSE,
              sanitize.text.function = function(x) paste0(x),
              hline.after = NULL, caption.placement = "top",
              add.to.row = list(pos = pos, command = labelsTable)
              )
    }

    if (transpose) {
        dotable(t(table), transpose)
    } else {
        if (ncols > 9) {
            ## This does not work for ncol = 9, 'Tuberculin' !!!
            ## One column in a table ?
            dotable(table[, 1:8])
            if (ncols > 16) {
                dotable(table[,  9:16])
                dotable(table[, 17:ncols])
            } else
                dotable(table[,  9:ncols])
        } else
            dotable(table)
    }
}
