`permalgorithm` <-
function (numSubjects, maxTime, Xmat, XmatNames = NULL, eventRandom = NULL, 
    censorRandom = NULL, betas, groupByD = FALSE) 
{
   if (length(betas) ==1 ) {Xmat = matrix(Xmat, ncol=1)}
   if ( dim(Xmat)[1] != numSubjects*maxTime)
       stop("length of Xmat does not equal numSubjects*maxTime")
   nc <- dim(Xmat)[2]   #The number of covariates
   I <- rep(seq(numSubjects), each = maxTime, times = nc)
   J <- rep(seq(maxTime), times = nc * numSubjects)
   K <- rep(seq(nc), each = maxTime * numSubjects)
   covArray <- array(0, dim = c(numSubjects, maxTime, nc))
   covArray[cbind(I, J, K)] <- Xmat

   if (is(eventRandom, "NULL")) 
     eventRandom <- function(n) sample(maxTime, numSubjects, replace = TRUE)
    if (is(eventRandom, "function")) {
        survivalTime <- as.integer(eventRandom(numSubjects))
    }
    else if (is(eventRandom, "numeric")) {
        if (length(eventRandom) != numSubjects) 
            stop("length of eventRandom is not equal to length of first dimension of Xmat")
        survivalTime <- as.integer(eventRandom)
    }
    else stop("eventRandom is neither numeric nor function")
    if (is(censorRandom, "NULL")) 
        censorRandom <- function(n) sample(maxTime, n, replace = TRUE)
    if (is(censorRandom, "function")) {
        censorTime <- as.integer(censorRandom(numSubjects))
    }
    else if (is(censorRandom, "numeric")) {
        if (length(censorRandom) != numSubjects) 
            stop("length of censorRandom is not equal to length of first dimension of Xmat")
        censorTime <- as.integer(censorRandom)
    }
    else stop("censorRandom is neither numeric nor function")
    if (min(survivalTime) <= 0) 
        stop("Not all event times are positive")
    if (!is(covArray, "array") || length(dim(covArray)) != 3) 
        stop("covArray is not a 3-dimensional array")
    if (length(betas) != dim(covArray)[3]) 
        stop("length of betas is not equal to length of third dimension of Xmat")
    notCensored <- ifelse(survivalTime <= apply(cbind(censorTime, 
        maxTime), MARGIN = 1, FUN = min), 1, 0)
    observedTime <- apply(cbind(survivalTime, censorTime, maxTime), 
        MARGIN = 1, FUN = min)
    I <- count <- integer(0)
    if (groupByD) {
        h <- 2 * observedTime + notCensored
        hBins = tabulate(h)
        I <- count <- integer(length(hBins))
        I[h] <- seq(along.with = h)
        I <- I[I > 0]
        count[I] <- hBins[h[I]]
    }
    else {
        I <- order(observedTime)
        count = rep(1, times = length(I))
    }
    p = .PermuteCovariateVectors(observedTime, notCensored, count, I, covArray, betas) 
    if (groupByD) {
        p[order(order(observedTime), sample(numSubjects))]
        J = order(h)
        tuples = cbind(obs.t = observedTime[J], d = notCensored[J], 
            id.tuples = J)
    }
    else {
        tuples = cbind(obs.t = observedTime[I], d = notCensored[I], 
            id.tuples = I)
    }
    info = data.frame(cbind(tuples, cov.id = p))
    ordered.info = info[order(info$cov.id), ]
    return(.form.data(ordered.info, maxTime, numSubjects, Xmat, 
        XmatNames))
}

