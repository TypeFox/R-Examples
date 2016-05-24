siarmcmcdirichletv4 <-
function (data, sources, corrections = 0, concdep = 0, iterations = 2e+05,
    burnin = 50000, howmany = 10000, thinby = 15, prior = rep(1,
        nrow(sources)), siardata = list(SHOULDRUN = FALSE))
{
    if (siardata$SHOULDRUN == FALSE) {
        siardata <- list()
        siardata$iterations <- iterations
        siardata$burnin <- burnin
        siardata$howmany <- howmany
        siardata$thinby <- thinby
        siardata$TITLE <- "SIAR data"
    }
    numsources <- nrow(sources)
    numdata <- nrow(data)
    numiso <- (ncol(sources) - 1)/2
    if (ncol(data) == numiso + 1) {
        data2 <- data.matrix(data[, 2:(numiso + 1)])   # # AJ add data.matrix and 2:(numiso + 1)
        numgroups <- max(data[, 1])
        startgroup <- as.vector(c(0, cumsum(table(data[, 1]))) +
            1)[1:numgroups]
        endgroup <- as.vector(cumsum(table(data[, 1])))
    }
    else {
        numgroups <- 1
        data2 <- data
        startgroup <- 1
        endgroup <- numdata
    }
    sourcenames <- as.character(sources[, 1])
    sourcedata <- sources[, 2:(2 * numiso + 1)]
    parameters <- matrix(1, ncol = (numsources + numiso) * numgroups,
        nrow = (siardata$iterations - siardata$burnin)/siardata$thinby)
    if (!is.data.frame(corrections)) {
        correctionsdata <- matrix(0, ncol = 2 * numiso, nrow = numsources)
    }
    else {
        correctionsdata <- matrix(as.double(as.matrix(corrections[,
            2:(2 * numiso + 1)])), nrow = numsources)
    }
    if (!is.data.frame(concdep)) {
        concdepdata <- matrix(1, ncol = 2 * numiso, nrow = numsources)
    }
    else {
        concdepdata <- matrix(as.double(as.matrix(concdep[, 2:(2 *
            numiso + 1)])), nrow = numsources)
    }
    BAD <- FALSE
    if (round((siardata$iterations - siardata$burnin)/siardata$thinby) !=
        (siardata$iterations - siardata$burnin)/siardata$thinby) {
        cat("Error in iterations, burnin or thinby: (iterations-burnin)/thinby must be an integer. \n \n")
        BAD <- TRUE
    }
    if (!is.numeric(data2[, 1])) {
        cat("Error in the target file - check this is numeric. \n \n")
        BAD <- TRUE
    }
    if (!is.numeric(sourcedata[, 1])) {
        cat("Error in the sources file - check this is numeric. \n \n")
        BAD <- TRUE
    }
    if (!is.numeric(correctionsdata[, 1])) {
        cat("Error in the corrections file - check this is numeric. \n \n")
        BAD <- TRUE
    }
    if (!is.numeric(concdepdata[, 1])) {
        cat("Error in the concentration dependence file - check this is numeric. \n \n")
        BAD <- TRUE
    }
    if (numgroups > 32) {
        cat("Error: too many groups - take some of the groups out before running. \n")
        BAD <- TRUE
    }
        
    if (BAD == FALSE) {
#         tempout <- .C("siarmcmcv4", as.integer(numdata), as.integer(numsources),
#             as.integer(numiso), as.integer(numgroups), as.integer(startgroup),
#             as.integer(endgroup), as.integer(siardata$iterations),
#             as.integer(siardata$burnin), as.integer(siardata$howmany),
#             as.integer(siardata$thinby), as.double(prior), as.data.frame(data2),
#             as.data.frame(concdepdata), as.data.frame(sourcedata),
#             as.data.frame(correctionsdata), as.data.frame(parameters))
        tempout <- .C("siarmcmcv4", as.integer(numdata), as.integer(numsources),
                      as.integer(numiso), as.integer(numgroups), as.integer(startgroup),
                      as.integer(endgroup), as.integer(siardata$iterations),
                      as.integer(siardata$burnin), as.integer(siardata$howmany),
                      as.integer(siardata$thinby), as.double(prior), as.double(data2),
                      as.double(concdepdata), as.double(as.matrix(sourcedata)),
                      as.double(correctionsdata), as.double(parameters))
        parameters <- matrix(tempout[[16]], ncol = (numsources + numiso) * numgroups,
                             nrow = (siardata$iterations - siardata$burnin)/siardata$thinby)
        if (numgroups == 1) {
            colnames(parameters) <- c(sourcenames, paste("SD",
                seq(1, numiso), sep = ""))
        }
        else {
            colnames(parameters) <- paste(rep(paste(c(sourcenames,
                paste("SD", seq(1, numiso), sep = "")), "G",
                sep = ""), times = numgroups), sort(rep(seq(1,
                numgroups), times = numsources + numiso)), sep = "")
        }
        return(list(targets = data, sources = sources, corrections = corrections,
            concdep = concdep, PATH = siardata$PATH, TITLE = siardata$TITLE,
            numgroups = tempout[[4]], numdata = tempout[[1]],
            numsources = tempout[[2]], numiso = tempout[[3]],
            SHOULDRUN = TRUE, GRAPHSONLY = FALSE, EXIT = FALSE,
            SIARSOLO = FALSE, output = parameters))
    }
    else {
        cat("Problems with inputs: siar has not been run. \n")
        return(siardata)
    }
}
