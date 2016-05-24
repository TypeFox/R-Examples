getConTable <- function(actual, pre){
    if(length(actual) != length(pre))
        stop("The length of 'actual' is not equal to the length of 'pre'.")
    aclass <- unique(actual)
    pclass <- unique(pre)
    if(! (all(pclass %in% aclass)))
        stop("The range of the predicted classes is out of the range of the actual classes.")
    class <- sort(aclass)
    nclass <- length(class)
    table <- matrix(0, nrow=nclass, ncol=nclass)
    for(i in 1:nclass){
        focus <- pre[actual == i]
        n <- length(focus)
        for(j in 1:nclass){
            table[i,j] <- sum(focus==j)/n
        }
    }
    t(table)
}

getDSM <- function(actual, pre){
    if(length(actual) != length(pre))
        stop("The length of 'actual' is not equal to the length of 'pre'.")
    aclass <- unique(actual)
    pclass <- unique(pre)
    if(! (all(pclass %in% aclass)))
        stop("The range of the predicted classes is out of the range of the actual classes.")
    class <- sort(aclass)
    nclass <- length(class)
    DSM <- rep(0, nclass)
    for(i in 1:nclass){
        comp <- pre == i
        true <- actual == i
        inter <- sum(comp & true)
        DSM[i] <- 2 * inter / (sum(comp) + sum(true))
    }
    DSM
}


measureMRI <- function(intvec=NULL, actual, pre){
    if (!all(dim(actual) == dim(pre)))
        stop("The dimension of 'actual' does not match that of the 'pre'.")
    if (!is.null(intvec) && !is.vector(intvec))
        stop("'intvec' has to be a vector.")
    if (!is.null(intvec) && length(intvec) != nrow(pre))
        stop("The number of intensity values does not match the dimension of the 'pre'.")
    actual.discrete <- max.col(actual)
    pre.discrete <- max.col(pre)

    #aclass <- unique(actual.discrete)
    #pclass <- unique(pre.discrete)
    #if (! (all(pclass %in% aclass)))
    #    stop("The range of the predicted classes is out of the range of the actual classes.")

    mse <- mean((pre - actual) ^ 2)
  
    pvolume <- unlist(lapply(1:ncol(actual), function(i) sum(pre.discrete==i)))
    avolume <- unlist(lapply(1:ncol(actual), function(i) sum(actual.discrete==i)))
    rseVolume <-  abs(pvolume - avolume) / avolume
                   
    misclass <- mean(actual.discrete != pre.discrete)

    DSM <- getDSM(actual.discrete, pre.discrete)

    conTable <- getConTable(actual.discrete, pre.discrete)

    if(! is.null(intvec)){
        intensity <- rep(intvec, 2)
        class <- c(actual.discrete, pre.discrete)
        g <- rep(c("actual", "predicted"), each=length(intensity)/2)
        dp <- densityplot(~ intensity | factor(class), groups = g,
                    plot.points = FALSE, ref = TRUE,
                    auto.key = list(columns = 2),
                    layout=c(1,ncol(pre)))
        plot(dp)
    }

    list(mse=mse, misclass=misclass, rseVolume=rseVolume, DSM=DSM, conTable=conTable)
}
