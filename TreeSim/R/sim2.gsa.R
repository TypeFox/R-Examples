sim2.gsa<-function (treelist, n, sampling) 
{
    treearray <- list()
    numbsvector <- list()
    timesvector <- list()
    timeperiodvector <- list()
    timeoverallvector <- vector()
    sumtimes <- 0
    k <- 0
    for (j in 1:length(treelist)) {
        if (class(treelist[[j]]) == "phylo") {
            numbs <- getnumbs(treelist[[j]])
            temp <- gettimelength(numbs, n)
            numbsvector <- c(numbsvector, list(numbs))
            timesvector <- c(timesvector, list(temp))
            timeperiodvector <- c(timeperiodvector, list(gettimeperiod(numbs, 
                n)))
            timeoverallvector <- c(timeoverallvector, sum(temp))
            sumtimes <- sumtimes + sum(temp)
            k <- k + 1
        }
    }
    timemax <- sumtimes/k * sampling
    j <- 0
    for (k in 1:length(treelist)) {
        if (class(treelist[[k]]) == "phylo") {
            j <- j + 1
            timeperiods <- timeperiodvector[[j]]
            if (length(timeperiods) > 0) {
                tree <- treelist[[j]]
                timelength <- timesvector[[j]]
                timeoverall <- timeoverallvector[[j]]
                numbs <- numbsvector[[j]]
                rsamp <- runif(1, min = 0, max = 1)
                if (rsamp <= (timeoverall/timemax)) {
                  r <- runif(1, min = 0, max = timeoverall)
                  timepassed <- 0
                  jtest <- 0
                  while (r > timepassed) {
                    jtest <- jtest + 1
                    timepassed <- timepassed + timelength[jtest]
                  }
                  cutinterval <- timeperiods[jtest]
                  cuttime <- numbs[cutinterval, 2] - runif(1, 
                    min = 0, max = timelength[jtest])
                  if (cuttime > max(branching.times.complete(tree))) {
                    stop("error")
                  }
                  treecut <- cuttree(tree, cuttime)
                  treearray <- c(treearray, list(reorder(treecut)))
                }
            }
        }
    }
    treearray
}


