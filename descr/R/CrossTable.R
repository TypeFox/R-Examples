# These functions were developed from the function CrossTable of the package
# gmodels.  The original function had the following comments:
#
# Revision 2.2 2006/05/02
# Fix a bug when a matrix is passed as the 'x' argument
# Reported by Prof. Albert Sorribas same day
# Fix involved creating default values for RowData and ColData
# when there are no dimnames for the matrix

# Revision 2.1 2005/06/26
# Added 'dnn' argument to enable specification of dimnames
# as per table()
# Correct bug in SPSS output for 1d table, where proportions
# were being printed and not percentages ('%' output)

# Revision 2.0 2005/04/27
# Added 'format = "d"' to all table count output
# so that large integers do not print in
# scientific notation

GetDigitsList <- function(x)
{
    if(is.list(x$digits))
        dgts <- list(expected = ifelse(is.null(x$digits$expected), 1, x$digits$expected),
                     prop = ifelse(is.null(x$digits$prop), 1, x$digits$prop),
                     percent = ifelse(is.null(x$digits$percent), 1, x$digits$percent),
                     others = ifelse(is.null(x$digits$others), 1, x$digits$others))
    else
        dgts <- list(expected = x$digits[1], prop = x$digits[1],
                     percent = x$digits[1],  others = x$digits[1])
    dgts
}

CreateNewTab <- function(x, ...)
{
    nr <- dim(x$tab)[1]
    nc <- dim(x$tab)[2]
    nt <- cbind(x$tab, x$rs)
    colnames(nt)[ncol(nt)] <- gettext("Total", domain = "R-descr")

    dgts <- GetDigitsList(x)

    if(x$format == "SPSS"){
        hdd <- 100
        dgts$prop <- dgts$percent
    } else {
        hdd <- 1
        x$percent <- FALSE
    }

    appendlines <- function(nt, xx, rowlab, prct = FALSE, hasttl = FALSE)
    {
        if(prct)
            for(i in 1:nrow(xx))
                for(j in 1:ncol(xx))
                    xx[i, j] <- paste0(xx[i, j], "%")
        if(!hasttl)
            xx <- cbind(xx, rep("", nr))
        if(!x$row.labels)
            rowlab <- " "
        rownames(xx) <- rep(rowlab, nrow(xx))
        n <- dim(nt)[1] / nr
        nt <- rbind(nt, xx)
        idx <- integer()
        k <- 1
        l <- nr * n + 1
        for(i in 1:nr){
            for(j in 1:n){
                idx <- c(idx, k)
                k <- k + 1
            }
            idx <- c(idx, l)
            l <- l + 1
        }
        nt <- nt[idx, ]
        nt
    }

    if(x$expected){
        xx <- x$CST$expected
        xx <- format(round(xx, dgts$expected), trim = TRUE, ...)
        nt <- appendlines(nt, xx, gettext("expected", domain = "R-descr"))
    }

    if(x$prop.chisq){
        xx <- ((x$CST$expected - x$tab) ^ 2) / x$CST$expected
        xx <- format(round(xx, digits = dgts$others), trim = TRUE, ...)
        nt <- appendlines(nt, xx, gettext("chisq", domain = "R-descr"))
    }

    if(!is.na(x$prop.row[1])){
        xx <- cbind(x$prop.row, x$rs / x$gt)
        xx <- format(round(xx * hdd, digits = dgts$prop), trim = TRUE, ...)
        if(hdd == 100)
            nt <- appendlines(nt, xx, gettext("row %", domain = "R-descr"),
                              x$percent, TRUE)
        else
            nt <- appendlines(nt, xx, gettext("row prop.", domain = "R-descr"),
                              x$percent, TRUE)
    }

    if(!is.na(x$prop.col[1])){
        xx <- format(round(x$prop.col * hdd, digits = dgts$prop), trim = TRUE, ...)
        if(hdd == 100)
            nt <- appendlines(nt, xx, gettext("col %", domain = "R-descr"),
                              x$percent)
        else
            nt <- appendlines(nt, xx, gettext("col prop.", domain = "R-descr"),
                              x$percent)
    }

    if(!is.na(x$prop.tbl[1])){
        xx <- format(round(x$prop.tbl * hdd, digits = dgts$prop), trim = TRUE, ...)
        if(hdd == 100)
            nt <- appendlines(nt, xx, gettext("table %", domain = "R-descr"),
                              x$percent)
        else
            nt <- appendlines(nt, xx, gettext("table prop.", domain = "R-descr"),
                              x$percent)
    }

    if(!is.na(x$resid) && x$resid == TRUE){
        xx <- x$tab - x$CST$expected
        xx <- format(round(xx, digits = dgts$others), trim = TRUE, ...)
        nt <- appendlines(nt, xx, gettext("residual", domain = "R-descr"))
    }

    if(!is.na(x$sresid) && x$sresid == TRUE){
        xx <- x$CST$residual
        xx <- format(round(xx, digits = dgts$others), trim = TRUE, ...)
        nt <- appendlines(nt, xx, gettext("std. res.", domain = "R-descr"))
    }

    if(!is.na(x$asr[1])){
        xx <- format(round(x$asr, digits = dgts$others), trim = TRUE, ...)
        nt <- appendlines(nt, xx, gettext("adj. std. res.", domain = "R-descr"))
    }

    if(x$total.c){
        nt <- rbind(nt, c(x$cs, x$gt))
        rownames(nt)[nrow(nt)] <- gettext("Total", domain = "R-descr")

        # Add final row if necessary
        if(!is.na(x$prop.col[1])){
            xx <- format(round(hdd * x$cs / x$gt, digits = dgts$prop), trim = TRUE, ...)
            if(hdd == 100 && x$percent)
                xx <- paste0(xx, "%")
            nt <- rbind(nt, c(xx, ""))
        }
    }

    # Delete last column if necessary
    if(!x$total.r)
        nt <- nt[, 1:(ncol(nt)-1)]

    tdim <- dim(nt)
    tdnn <- list(rownames(nt), colnames(nt))
    names(tdnn) <- c(attr(nt, "RowData"), attr(nt, "ColData"))
    attributes(nt) <- NULL
    attr(nt, "dim") <- tdim
    dimnames(nt) <- tdnn
    nt
}

Abbrev1 <- function(x, lmt)
{
    if(nchar(x) <= lmt)
        return(x)
    xc <- charToRaw(x)
    len1 <- length(xc)

    # Two consecutive spaces
    i <- 1
    while(i < len1){
        if(xc[i] == 0x20 && xc[i+1] == 0x20)
            return(rawToChar(xc[-i]))
        i <- i + 1
    }

    # Trailing space or tab
    if(xc[len1] == 0x20 || xc[len1] == 0x09)
        return(rawToChar(xc[-len1]))

    # Last lower case ascii vowel (unless it is the first letter in a word)
    i <- len1
    while(i > 1){
        if(xc[i-1] != 0x20 && (xc[i] == 0x61 || xc[i] == 0x65 || xc[i] == 0x69 || xc[i] == 0x6f || xc[i] == 0x75)){
            xc <- xc[-i]
            return(rawToChar(xc))
        }
        i <- i - 1
    }

    # Last lower case ascii letter (unless it is the first letter in a word)
    i <- len1
    while(i > 0){
        if(xc[i-1] != 0x20 && xc[i] > 0x60 && xc[i] < 0x7b){
            xc <- xc[-i]
            return(rawToChar(xc))
        }
        i <- i - 1
    }

    # Last letter
    x <- unlist(strsplit(x, ""))
    len <- length(x)
    x <- x[-len]
    paste0(x, collapse = "")
}

CrossTable <- function (x, y,
    digits = list(expected = 1, prop = 3, percent = 1, others = 3),
    max.width = NA, expected = FALSE,
    prop.r = TRUE, prop.c = TRUE, prop.t = TRUE, prop.chisq = TRUE,
    chisq = FALSE, fisher = FALSE, mcnemar = FALSE, resid = FALSE,
    sresid = FALSE, asresid = FALSE, missing.include = FALSE,
    drop.levels = TRUE, format = c("SAS", "SPSS"), dnn = NULL,
    cell.layout = TRUE, row.labels = !cell.layout,
    percent = (format == "SPSS" && !row.labels),
    total.r, total.c, xlab = NULL, ylab = NULL, ...)
{
    format = match.arg(format)

    RowData <- deparse(substitute(x))
    if (missing(y))
        ColData <- NA
    else
        ColData <- deparse(substitute(y))

    ## Ensure that max.width >= 1
    if (!is.na(max.width) && max.width < 1)
        stop("max.width must be >= 1")
    ## Set 'x' vector flag
    vector.x <- FALSE

    if (missing(y))
    {
        ## is x a vector?
        if (is.null(dim(x)))
        {
            if (missing.include)
                x <- no.drop.levels(x)
            if (drop.levels)
                x <- factor(x)
            tab <- t(as.matrix(table(x)))
            vector.x <- TRUE
        }
        ## is x a matrix?
        else if (length(dim(x) == 2))
        {
            if(any(x < 0) || any(is.na(x)))
                stop("all entries of x must be nonnegative and finite")

            ## Check to see if x has names(dimnames) defined. If yes, use these for
            ## 'RowData' and 'ColData' labels, else create blank ones
            ## This can be overridden by setting 'dnn' values
            if (is.null(names(dimnames(x))))
            {
                RowData <- ""
                ColData <- ""
            } else {
                RowData <- names(dimnames(x))[1]
                ColData <- names(dimnames(x))[2]
            }

            ## Add generic column and rownames if required
            ## check each separately, in case user has defined one or the other
            if (is.null(rownames(x)))
                rownames(x) <- paste("[", 1:nrow(x), ",]", sep = "")
            if (is.null(colnames(x)))
                colnames(x) <- paste("[,", 1:ncol(x), "]", sep = "")

            tab <- x
        }
        else
            stop("x must be either a vector or a 2 dimensional matrix, if y is not given")
    } else {
        if(missing.include){
            x <- no.drop.levels(x)
            y <- no.drop.levels(y)
        } 
        if(drop.levels){
            x <- factor(x)
            y <- factor(y)
        }

        if(length(x) != length(y))
            stop("x and y must have the same length")

        ## Generate table
        tab <- table(x, y)
    }

    ## Create Titles for Table From Vector Names
    ## At least 2 x 2 table only (for now)
    if (all(dim(tab) >= 2)) {
        if (!is.null(dnn)) {
            if (length(dnn) != 2)
                stop("dnn must have length of 2, one element for each table dimension")
            else {
                RowData <- dnn[1]
                ColData <- dnn[2]
            }
        }
    }

    ## if tab is not at least a 2 x 2, do not do stats
    ## even if any set to TRUE. Do not do col/table props
    if (any(dim(tab) < 2))
        prop.c <- prop.chisq <- chisq <- expected <- fisher <- mcnemar <- FALSE

    if (vector.x && dim(tab)[2] < 2)
        prop.r <- FALSE
    if (!vector.x && dim(tab)[1] < 2)
        prop.r <- FALSE

    CPR <- CPC <- CPT <- GT <- RS <- TotalN <- CSTc <- CST <- ASR <- FTt <-
        FTl <- FTg <- McN <- McNc <- NA

    ## Generate cell proportion of row
    if(prop.r)
        CPR <- prop.table(tab, 1)

    ## Generate cell proportion of col
    if(prop.c)
        CPC <- prop.table(tab, 2)

    ## Generate cell proportion of total
    if(prop.t)
        CPT <- prop.table(tab)

    ## Generate summary counts
    GT <- sum(tab)
    RS <- rowSums(tab)
    CS <- colSums(tab)

    if (length(dim(x) == 2))
        TotalN <- GT
    else
        TotalN <- length(x)

    ## Perform Chi-Square Tests
    if (expected || chisq || prop.chisq || resid || sresid || asresid) {
        if(!chisq && !prop.chisq){
            wv <- getOption("warn")
            options(warn = -1)
        }
        CST <- chisq.test(tab, correct = FALSE, ...)
        if (all(dim(tab) == 2))
            CSTc <- chisq.test(tab, correct = TRUE, ...)
        if(!chisq && !prop.chisq){
            options(warn = wv)
        }
    }

    if (asresid & !vector.x)
        ASR <- (CST$observed-CST$expected)/sqrt(CST$expected*((1-RS/GT) %*% t(1-CS/GT)))

    if (fisher)
    {
        try(FTt <- fisher.test(tab, alternative = "two.sided"))
        if (all(dim(tab) == 2))
        {
            FTl <- fisher.test(tab, alternative = "less")
            FTg <- fisher.test(tab, alternative = "greater")
        }
    }

    if (mcnemar)
    {
        if(dim(tab)[1] == dim(tab)[2])
            McN <- mcnemar.test(tab, correct = FALSE)
        if (all(dim(tab) == 2))
            McNc <- mcnemar.test(tab, correct = TRUE)
    }

    if(!missing(total.r)){
        if(!is.logical(total.r))
            stop(gettext("total.r must be logical", domain = "R-descr"))
        if(missing(total.c))
            total.c <- total.r
    }
    if(!missing(total.c)){
        if(!is.logical(total.c))
            stop(gettext("total.c must be logical", domain = "R-descr"))
        if(missing(total.r))
            total.r <- total.c
    }
    if(missing(total.r) & missing(total.c))
        total.r <- total.c <- TRUE

    res <- list(tab = tab, prop.row = CPR, prop.col = CPC, prop.tbl = CPT,
                gt = GT, rs = RS, cs = CS, total.n = TotalN, chisq = chisq,
                CST = CST, chisq.corr = CSTc, fisher.ts = FTt,
                fisher.lt = FTl, fisher.gt = FTg, print.mcnemar = mcnemar,
                mcnemar = McN, mcnemar.corr = McNc, asr = ASR,
                RowData = RowData, ColData = ColData, digits = digits,
                max.width = max.width, vector.x = vector.x,
                expected = expected, prop.chisq = prop.chisq, resid = resid,
                sresid = sresid, asresid = asresid, format = format,
                cell.layout = cell.layout, row.labels = row.labels,
                percent = percent, total.r = total.r, total.c = total.c)

    # Add "t" element to avoid breaking pander package which manipulates the
    # CrossTable object:
    res$t <- res$tab

    # Attributes for plotting
    attr(res, "xlab") <- xlab
    attr(res, "ylab") <- ylab

    class(res) <- "CrossTable"
    res
}

print.CrossTable <- function(x, ...)
{
    argl <- list(...)
    for(n in names(argl))
        if(n %in% names(x))
            x[[n]] <- argl[[n]]

    nt <- CreateNewTab(x, ...)
    tab <- x$tab
    CPR <- x$prop.row
    CPC <- x$prop.col
    CPT <- x$prop.tbl
    GT <- x$gt
    RS <- x$rs
    CS <- x$cs
    TotalN <- x$total.n
    chisq <- x$chisq
    CST <- x$CST
    CSTc <- x$chisq.corr
    FTt <- x$fisher.ts
    FTl <- x$fisher.lt
    FTg <- x$fisher.gt
    McN <- x$mcnemar
    McNc <- x$mcnemar.corr
    ASR <- x$asr
    RowData <- x$RowData
    ColData <- x$ColData
    max.width <- x$max.width
    vector.x <- x$vector.x
    expected <- x$expected
    prop.r <- !is.na(CPR[1])
    prop.c <- !is.na(CPC[1])
    prop.t <- !is.na(CPT[1])
    prop.chisq <- x$prop.chisq
    fisher <- (class(FTt) == "htest")
    resid <- x$resid
    sresid <- x$sresid
    asresid <- x$asresid
    mcnemar <- x$print.mcnemar
    format <- x$format
    cell.layout <- x$cell.layout
    total.r <- x$total.r
    total.c <- x$total.c
    outDec <- getOption("OutDec")

    dgts <- GetDigitsList(x)

    if(format == "SAS") {
        hdd <- 1
    } else {
        if (format == "SPSS") {
            hdd <- 100
            dgts$prop <- dgts$percent
        } else {
            stop("unknown format")
        }
    }

    if(vector.x)
        expected <- prop.chisq <- prop.c <- prop.t <- resid <- sresid <- asresid <- FALSE

    #### Printing the tables

    ## Print Cell Layout
    if(cell.layout){
        cat("  ", gettext("Cell Contents", domain = "R-descr"), "\n")
        cat("|-------------------------|\n")
        if (format=="SAS") {
            cat(gettext("|                       N |", domain = "R-descr"), "\n")
            if (expected)
                cat(gettext("|              Expected N |", domain = "R-descr"), "\n")
            if (prop.chisq)                                                     
                cat(gettext("| Chi-square contribution |", domain = "R-descr"), "\n")
            if (prop.r)                                                         
                cat(gettext("|           N / Row Total |", domain = "R-descr"), "\n")
            if (prop.c)                                                         
                cat(gettext("|           N / Col Total |", domain = "R-descr"), "\n")
            if (prop.t)                                                         
                cat(gettext("|         N / Table Total |", domain = "R-descr"), "\n")
        } else if (format == "SPSS") {
            cat(gettext("|                   Count |", domain = "R-descr"), "\n")
            if (expected)
                cat(gettext("|         Expected Values |", domain = "R-descr"), "\n")
            if (prop.chisq)                                                       
                cat(gettext("| Chi-square contribution |", domain = "R-descr"), "\n")
            if (prop.r)                                                           
                cat(gettext("|             Row Percent |", domain = "R-descr"), "\n")
            if (prop.c)                                                           
                cat(gettext("|          Column Percent |", domain = "R-descr"), "\n")
            if (prop.t)                                                           
                cat(gettext("|           Total Percent |", domain = "R-descr"), "\n")
        }
        if (resid)                                                            
            cat(gettext("|                Residual |", domain = "R-descr"), "\n")
        if (sresid)                                                           
            cat(gettext("|            Std Residual |", domain = "R-descr"), "\n")
        if (asresid)                                                          
            cat(gettext("|           Adj Std Resid |", domain = "R-descr"), "\n")
        cat("|-------------------------|\n")
    }

    ## Print 1 X N vector
    if (vector.x) {
        ## Set consistent column widths based upon dimnames and table values
        strt <- formatC(unclass(tab), digits = dgts$prop, format = "f", width = 0, decimal.mark = outDec)
        CWidth <- max(dgts$prop + 2, c(nchar(strt, type = "width"),
                                    nchar(dimnames(tab)[[2]], type = "width")))
        if(prop.r){
            if(vector.x)
                strt <- formatC(unclass(CPT), digits = dgts$prop, format = "f", width = 0, decimal.mark = outDec)
            CWidth <- max(CWidth, nchar(strt, type = "width"))
        }

        ## Create row separators
        RowSep <- paste(rep("-", CWidth + 2), collapse = "")

        if(is.na(max.width))
            max.width = floor((getOption("width") - 2) / (CWidth + 3))
        if (length(tab) > max.width) {
            ## set breakpoints for output based upon max.width
            final.row <- length(tab) %% max.width
            max <- length(tab) - final.row
            ## Define breakpoint indices for each row
            start <- seq(1, max, max.width)
            end <- start + (max.width - 1)
            ## Add final.row if required
            if (final.row > 0) {
                start <- c(start, end[length(end)] + 1)
                end <- c(end, end[length(end)] + final.row)
            }
        } else {
            ## Each value printed horizontally in a single row
            start <- 1
            end <- length(tab)
        }

        cat("\n")

        for (i in 1:length(start)) {
            cat("| ")
            cat(paste(formatC(dimnames(tab)[[2]][start[i]:end[i]], width = CWidth, format = "s"), collapse = " | "))
            cat(" |\n|")
            cat(rep(RowSep, (end[i] - start[i]) + 1), sep = "|")
            cat("|\n| ")
            cat(formatC(tab[, start[i]:end[i]], width = CWidth, format = "d"), sep = " | ")
            if(prop.r){
                cat(" |\n| ")
                cat(formatC(CPT[, start[i]:end[i]] * hdd, width = CWidth,
                            digits = dgts$prop, format = "f", decimal.mark = outDec), sep = " | ")
            }
            cat(" |\n|")
            cat(rep(RowSep, (end[i] - start[i]) + 1), sep = "|")
            cat("|\n\n")

        }  ## End of for (i in 1:length(start))

        if(format == "SPSS" && GT < TotalN)
            cat("\n", gettext("Number of Missing Observations:", domain = "R-descr"),
                " ", TotalN-GT, " (", 100*(TotalN-GT)/TotalN, "%)\n", sep = "")
        return(invisible(x))
    } ## End of if (vector.x)


    nr <- nrow(nt)
    nc <- ncol(nt)
    rnames <- rownames(nt)
    cnames <- colnames(nt)

    # Check column widths and fix them if necessary
    availablewidth <- getOption("width")
    minColWd <- max(nchar(nt, type = "width"))

    clabwidth <- nchar(cnames, type = "width")
    rlabwidth <- max(nchar(rnames, type = "width"), nchar(RowData, type = "width")) + 1
    totalwidth <- rlabwidth + (sum(clabwidth + 3))

    # Abbreviate row and column labels. Withdraw one char of the newlimit
    # label until the rows fit in the screen
    while(totalwidth > availablewidth){
        newlimit <- max(c(rlabwidth - 1, clabwidth)) - 1
        # FIXME: The algorithm should consider the minColWd of each column
        if(newlimit < minColWd)
            break

        rnames <- sapply(rnames, Abbrev1, newlimit)
        cnames <- sapply(cnames, Abbrev1, newlimit)
        RowData <- Abbrev1(RowData, newlimit)

        clabwidth <- nchar(cnames, type = "width")
        clabwidth[clabwidth < minColWd] <- minColWd
        rlabwidth <- max(nchar(rnames, type = "width"), nchar(RowData, type = "width")) + 1
        totalwidth <- rlabwidth + (sum(clabwidth + 3))
    }

    minColWd <- apply(nt, 2, function(x) max(nchar(x, type = "width")))
    clabwidth <- apply(cbind(clabwidth, minColWd), 1, max)
    dashedline <- rep("-", sum(clabwidth) + 3 * nc + rlabwidth)
    ddashedline <- gsub("-", "=", dashedline)

    # Calculate horizontal line locations
    if(x$total.c && !is.na(x$prop.col)[1])
        nrnt <- nrow(nt) - 2
    else if(x$total.c)
        nrnt <- nrow(nt) - 1
    else
        nrnt <- nrow(nt)
    n <- nrnt / nrow(tab)
    idxh <- seq(n+1, nrnt+1, n)
    # idxh <- idxh[idxh < nrow(nt)] # necessary when total.c = FALSE

    ## Print table cells
    cat("\n", ddashedline, "\n", sep = "")
    if(ColData != "")
        cat(formatC(" ", width = rlabwidth), "   ", ColData, "\n", sep = "", collapse = "")
    if(RowData == "")
        cat(formatC(" ", width = rlabwidth))
    else
        cat(formatC(RowData, width = rlabwidth, format = "s", flag = "-"))
    for(j in 1:nc)
        cat("  ", formatC(cnames[j], width = clabwidth[j]))
    cat("\n", dashedline, "\n", sep = "")
    for(i in 1:nr){
        if(i %in% idxh)
            cat(dashedline, "\n", sep = "")
        cat(formatC(rnames[i], width = rlabwidth, format = "s", flag = "-"), sep = "")
        if(x$percent && (prop.r || prop.c || prop.t)){
            for(j in 1:nc)
                if(grepl("%", nt[i, j]))
                    cat("   ", formatC(nt[i, j], width = clabwidth[j]), sep = "")
                else
                    cat("  ", formatC(nt[i, j], width = clabwidth[j]), " ", sep = "")
        } else {
            for(j in 1:nc)
                cat("  ", formatC(nt[i, j], width = clabwidth[j]))
        }
        cat("\n")
    }
    cat(ddashedline, "\n", sep = "")


    ## Print Statistics
    if (chisq)
    {
        cat("\n")
        cat(gettext("Statistics for All Table Factors", domain = "R-descr"),
            "\n\n", sep="")

        cat(CST$method, "\n")
        cat("------------------------------------------------------------\n")
        fp <- format.pval(CST$p.value, digits = dgts$others)
        pv <- paste("p", if(substr(fp, 1L, 1L) == "<") fp else paste("=", fp))
        cat(gettext("Chi^2 =", domain = "R-descr"), CST$statistic,
            "    ", gettext("d.f. =", domain = "R-descr"), CST$parameter,
            "    ", pv, "\n\n")

        if (all(dim(tab) == 2))
        {
            cat(CSTc$method, "\n")
            cat("------------------------------------------------------------\n")
            fp <- format.pval(CSTc$p.value, digits = dgts$others)
            pv <- paste("p", if(substr(fp, 1L, 1L) == "<") fp else paste("=", fp))
            cat(gettext("Chi^2 =", domain = "R-descr"), CSTc$statistic,
                "    ", gettext("d.f. =", domain = "R-descr"), CSTc$parameter,
                "    ", pv, "\n")
        }
    }

    ## Print McNemar tests
    if (is.na(McN[1]) == FALSE)
    {
        cat(rep("\n", 2))
        cat(McN$method, "\n")
        cat("------------------------------------------------------------\n")
        fp <- format.pval(McN$p.value, digits = dgts$others)
        pv <- paste("p", if(substr(fp, 1L, 1L) == "<") fp else paste("=", fp))
        cat(gettext("Chi^2 =", domain = "R-descr"), McN$statistic,
            "    ", gettext("d.f. =", domain = "R-descr"), McN$parameter,
            "    ", pv, "\n\n")

        if (is.na(McNc[1]) == FALSE)
        {
            cat(McNc$method, "\n")
            cat("------------------------------------------------------------\n")
            fp <- format.pval(McNc$p.value, digits = dgts$others)
            pv <- paste("p", if(substr(fp, 1L, 1L) == "<") fp else paste("=", fp))
            cat(gettext("Chi^2 =", domain = "R-descr"), McNc$statistic,
                "    ", gettext("d.f. =", domain = "R-descr"), McNc$parameter,
                "    ", pv, "\n")
        }
    }

    ## Pint Fisher Tests
    if (fisher)
    {
        cat(rep("\n", 2))

        cat(gettext("Fisher's Exact Test for Count Data", domain = "R-descr"))
        cat("\n------------------------------------------------------------\n")

        if (all(dim(tab) == 2))
        {
            cat(gettext("Sample estimate odds ratio:", domain = "R-descr"), FTt$estimate, "\n\n")

            cat(gettext("Alternative hypothesis: true odds ratio is not equal to 1",
                    domain = "R-descr"), "\n")
            fp <- format.pval(FTt$p.value, digits = dgts$others)
            pv <- paste("p", if(substr(fp, 1L, 1L) == "<") fp else paste("=", fp))
            cat(pv, "\n")
            cat(gettext("95% confidence interval:", domain = "R-descr"), FTt$conf.int, "\n\n")

            cat(gettext("Alternative hypothesis: true odds ratio is less than 1",
                    domain = "R-descr"), "\n")
            fp <- format.pval(FTl$p.value, digits = dgts$others)
            pv <- paste("p", if(substr(fp, 1L, 1L) == "<") fp else paste("=", fp))
            cat(pv, "\n")
            cat(gettext("95% confidence interval:", domain = "R-descr"), FTl$conf.int, "\n\n")

            cat(gettext("Alternative hypothesis: true odds ratio is greater than 1",
                    domain = "R-descr"), "\n")
            fp <- format.pval(FTg$p.value, digits = dgts$others)
            pv <- paste("p", if(substr(fp, 1L, 1L) == "<") fp else paste("=", fp))
            cat(pv, "\n")
            cat(gettext("95% confidence interval:", domain = "R-descr"), FTg$conf.int, "\n\n")
        }
        else
        {
            cat(gettext("Alternative hypothesis: two.sided", domain = "R-descr"), "\n")
            fp <- format.pval(FTt$p.value, digits = dgts$others)
            pv <- paste("p", if(substr(fp, 1L, 1L) == "<") fp else paste("=", fp))
            cat(pv, "\n")
        }
    } ## End Of If(Fisher) Loop

    #  cat(rep("\n", 2))

    if(format == "SPSS"){
        if (any(dim(tab) >= 2) & any(chisq, mcnemar, fisher))
        {
            MinExpF = min(CST$expected)
            cat("       ", gettext("Minimum expected frequency:", domain = "R-descr"), MinExpF, "\n")
            NMinExpF = length(CST$expected[which(CST$expected<5)])
            if (NMinExpF > 0)
            {
                NCells = length(CST$expected)
                cat(gettext("Cells with Expected Frequency < 5:", domain = "R-descr"),
                    " ", NMinExpF, " ", gettext("of", domain = "R-descr"), " ",
                    NCells, " (", 100*NMinExpF/NCells, "%)\n", sep = "")
            }
            cat("\n")

        } ## End of if (any(dim(tab)...))
    }
    return(invisible(x))
}

as.data.frame.CrossTable <- function(x, ...) as.data.frame(x$tab, ...)


# Needed by tableStyles() of odfWeave package:
dim.CrossTable <- function(x){
    dim(x$tab) + 1
}

