.invalidBlockDesign <- function(data,
                                response = "Response",
                                replicate = "Replicate",
                                sample = "Sample",
                                step = "Dilution",
                                sampleStep = "SampleStep"
                                ) {
    if (!any(dimnames(data)[[2]] == sampleStep)) {
        SampleStep <- paste0(as.character(unlist(data[sample])), ":",
                             as.character(unlist(data[step])))
        names(SampleStep) <- "SampleStep"
        data <- cbind(data, SampleStep = SampleStep)
    }
    dupletsRS <-  data[duplicated(data[c(replicate, sampleStep)]),]
    Error <- FALSE
    if (dim(dupletsRS)[1] > 0)
        Error <- TRUE
    tRS <- table(data[c(replicate, sampleStep)])
    Xreplicate <- NULL
    Xcolumn <- NULL
    XsampleStep <- NULL
    if (length(which(tRS != 1)) > 0) {
        Error <- TRUE
        Replicates  <- unique(sort(unlist(data[replicate])))
        SampleSteps <- unique(sort(unlist(data[sampleStep])))
        byReplicate <- split(data[, sampleStep], data[, replicate])
        bySample    <- split(data[, replicate],  data[, sampleStep])
        funReplicate      <- function(item)
            Replicates[is.na(match(Replicates, unlist(item)))]
        funSample   <- function(item)
            SampleSteps[is.na(match(SampleSteps, unlist(item)))]
        for(i in 1:length(byReplicate)) {
            item <- byReplicate[[i]]
            missing <- data.frame(sampleStep = funSample(item))
            if (dim(missing)[[1]] > 0)
                Xreplicate <- rbind(Xreplicate,
                                    cbind(replicate = rep(i,
                                              dim(missing)[1]), missing))
        }
        for(i in 1:length(bySample)) {
            item <- bySample[[i]]
            missing <- data.frame(replicate = funReplicate(item))
            if (dim(missing)[[1]] > 0)
                XsampleStep <- rbind(XsampleStep,
                                     cbind(sampleStep = rep(SampleSteps[i],
                                               dim(missing)[1]), missing))
        }
    }
    Result <- NULL
    if (Error)
        Result <- list(
            countsReplicateSample = tRS,
            dupletsReplicateSample = dupletsRS,
            missingByReplicate = Xreplicate,
            missingBySampleStep = XsampleStep)
    return(Result)
}
