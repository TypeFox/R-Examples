.invalidLatinSquare <- function(data,
                                response = "Response",
                                row = "Row",
                                column = "Column",
                                sample = "Sample",
                                step = "Dilution",
                                sampleStep = "SampleStep"
                                ) {
    ## mc <- function(sample, step)
    ##     paste0("", sample, "\\textsubscript{", step, "}")
    if (!any(dimnames(data)[[2]] == sampleStep)) {
        SampleStep <- paste0(as.character(unlist(data[sample])), ":",
                             as.character(unlist(data[step])))
        ## SampleStep <- mc(as.character(unlist(data[sample])),
        ##                 as.character(unlist(data[step])))
        names(SampleStep) <- "SampleStep"
        data <- cbind(data, SampleStep = SampleStep)
    }
    dupletsRC <-  data[duplicated(data[c(row, column)]),]
    dupletsRS <-  data[duplicated(data[c(row, sampleStep)]),]
    dupletsCS <-  data[duplicated(data[c(column, sampleStep)]),]
    Error <- FALSE
    if (dim(dupletsRC)[1] > 0)
        Error <- TRUE
    if (dim(dupletsRS)[1] > 0)
        Error <- TRUE
    if (dim(dupletsCS)[1] > 0)
        Error <- TRUE
    tRC <- table(data[c(row, column)])
    tRS <- table(data[c(row, sampleStep)])
    tCS <- table(data[c(column, sampleStep)])
    Xrow <- NULL
    Xcolumn <- NULL
    XsampleStep <- NULL
    if (length(which(tRC != 1)) > 0 |
        length(which(tRS != 1)) > 0 |
        length(which(tCS != 1)) > 0) {
        Error <- TRUE
        Rows        <- unique(sort(unlist(data[row])))
        Columns     <- unique(sort(unlist(data[column])))
        SampleSteps <- unique(sort(unlist(data[sampleStep])))
        byRow       <- split(data[, c(sampleStep, column)], data[, row])
        byColumn    <- split(data[, c(sampleStep, row)], data[, column])
        bySample    <- split(data[, c(row, column)], data[, sampleStep])
        funRow      <- function(item)
            Rows[is.na(match(Rows, unlist(item[row])))]
        funCol      <- function(item)
            Columns[is.na(match(Columns, unlist(item[column])))]
        funSample   <- function(item)
            SampleSteps[is.na(match(SampleSteps,
                                    unlist(item[sampleStep])))]
        for(i in 1:length(byRow)) {
            item <- byRow[[i]]
            missing <- data.frame(column = funCol(item),
                                  sampleStep = funSample(item))
            if (dim(missing)[[1]] > 0)
                Xrow <- rbind(Xrow,
                              cbind(row = rep(i, dim(missing)[1]), missing))
        }
        for(i in 1:length(byColumn)) {
            item <- byColumn[[i]]
            missing <- data.frame(row = funRow(item),
                                  sampleStep = funSample(item))
            if (dim(missing)[[1]] > 0)
                Xcolumn <- rbind(Xcolumn,
                                 cbind(column = rep(i,
                                           dim(missing)[1]), missing))
        }
        for(i in 1:length(bySample)) {
            item <- bySample[[i]]
            missing <- data.frame(row = funRow(item), column = funCol(item))
            if (dim(missing)[[1]] > 0)
                XsampleStep <- rbind(XsampleStep,
                                     cbind(sampleStep = rep(SampleSteps[i],
                                               dim(missing)[1]), missing))
        }
    }
    Result <- NULL
    if (Error)
        Result <- list(
            countsRowColumn = tRC,
            countsRowSample = tRS,
            countsColumnSample = tCS,
            dupletsRowColumn = dupletsRC,
            dupletsRowSample = dupletsRS,
            dupletsColumnSample = dupletsCS,
            missingByRow = Xrow,
            missingByColumn = Xcolumn,
            missingBySampleStep = XsampleStep)
    return(Result)
}
