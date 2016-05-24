searchrules <- function(d,
                        lhs=2:ncol(d),
                        rhs=1,
                        tnorm=c("goedel", "goguen", "lukasiewicz"),
                        n=100,
                        best=c("confidence"),
                        minSupport=0.02,
                        minConfidence=0.75,
                        maxConfidence=1,
                        maxLength=4,
                        numThreads=1,
                        trie=(maxConfidence < 1)) {

    if (!is.fsets(d)) {
        stop("'d' must be an instance of class 'fsets'")
    }
    if (ncol(d) < 2) {
        stop("'d' must have at least 2 columns")
    }
    if (nrow(d) < 1) {
        stop("'d' must not be empty")
    }

    if (!is.vector(lhs, mode='numeric')) {
        stop("'lhs' must be numeric vector")
    }
    if (min(lhs) < 1 || max(lhs) > ncol(d)) {
        stop("'lhs' must contain valid indexes of the columns of 'd'")
    }

    if (!is.vector(rhs, mode='numeric')) {
        stop("'rhs' must be numeric vector")
    }
    if (min(rhs) < 1 || max(rhs) > ncol(d)) {
        stop("'rhs' must contain valid indexes of the columns of 'd'")
    }

    if (!is.numeric(n) || length(n) > 1 || n < 0) {
        stop("'n' must be non-negative number (zero means unlimited number of rules)")
    }

    if (minSupport < 0 || minSupport > 1) {
        stop("'minSupport' must be value from interval [0, 1]")
    }

    if (minConfidence < 0 || minConfidence > 1) {
        stop("'minConfidence' must be value from interval [0, 1]")
    }
    
    if (maxConfidence < 0 || maxConfidence > 1) {
        stop("'maxConfidence' must be value from interval [0, 1]")
    }
    if (maxConfidence < minConfidence) {
        stop("'maxConfidence' must be greater or equal to 'minConfidence'")
    }
    
    if (maxLength < 0) {
        maxLength <- -1
    }

    if (!is.numeric(numThreads) || length(numThreads) > 1 || numThreads < 1) {
        stop("'numThreads' must be positive integer number")
    }

    tnorm <- match.arg(tnorm)
    if (tnorm == 'goedel') {
        tnorm <- 'minimum'
    } else if (tnorm == 'goguen') {
        tnorm <- 'product'
    }

    best = match.arg(best)

    config <- list(vars=as.numeric(as.factor(vars(d)[colnames(d)])),
                   minSupport=minSupport,
                   minConfidence=minConfidence,
                   maxConfidence=maxConfidence,
                   maxLength=maxLength,
                   lhs=lhs - 1,
                   rhs=rhs - 1,
                   tnorm=tnorm,
                   n=n,
                   best=best,
                   numThreads=numThreads,
                   trie=trie)

    #str(config)

    result <- .Call("search", d, config, PACKAGE="lfl")

    map <- colnames(d)
    names(map) <- 1:ncol(d)

    rules <- lapply(result$rules, function(x) {
                        x <- x + 1
                        x[] <- map[unlist(x)]
                        return(x)
                    })
    if (is.null(result$statistics) || length(result$statistics) <= 0) {
        stats <- matrix(0, nrow=0, ncol=0)
    } else {
        stats <-matrix(unlist(result$statistics),
                       byrow=TRUE,
                       nrow=length(result$statistics))
        colnames(stats) <- c('support', 'lhsSupport', 'rhsSupport', 'confidence', 'lift', 'loLift', 'hiLift')
        #volume <- .computeVolume(d, rules, tNorm(tnorm))
        #cat("Rules:\n")
        #print(rules)
        #cat("Stats:\n")
        #str(stats)
        #cat("volume:\n")
        #str(volume)
        #stats <- cbind(stats, volume=volume)
    }
    res <- farules(rules=rules,
                   statistics=stats)

    #if (maxConfidence < 1) {
        #cat('Post-filtering maxconf on rules: \n')
        #print(length(res$rules))
        #res <- maxconf(res, maxConfidence)
        #cat('After filter: ')
        #print(length(res$rules))
    #}
    return(res)
}

