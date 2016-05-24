


#this programming used unequal-variance t-test, if use equal-variance t-test, can get little different genelists.
# i.e. in function entity.build2, method = "t.equalvar" instead of "t")
# ref library(multtest)  function: mt.teststat
# If test="t", the tests are based on two-sample Welch t-statistics (unequal variances).
# If test="t.equalvar", the tests are based on two-sample t-statistics with equal variance for the two samples. The square of the t-statistic is equal to an F-statistic for k=2.



entitybuild2 <- function (expr.mat, ALLtype = NULL, type, dataset = NULL, minSampleNum = 3, 
    method = "t", random = FALSE) 
{
    n.patient <- length(ALLtype)
    #cat("entities building")
    if (ncol(expr.mat) != length(ALLtype)) 
        stop("expr.mat and ALLtype not match!")
    n.type <- length(type)
    entity.labels <- array(0, dim = ncol(expr.mat))
    for (i in 1:n.type) entity.labels[which(ALLtype == type[i])] <- i
    if (n.type == 2) {
        i <- 1
        j <- 2
        if (length(which(entity.labels == i)) > minSampleNum) 
            if (length(which(entity.labels == j)) > minSampleNum) {
                #cat(".")
                entity1 <- paste(dataset, type[i], sep = "")
                entity2 <- paste(dataset, type[j], sep = "")
                groupA <- which(entity.labels == i)
                nA <- length(groupA)
                if (random) {
                  groupA <- sample(n.patient, nA)
                }
                groupB <- which(entity.labels == j)
                nB <- length(groupB)
                if (random) {
                  groupB <- (1:n.patient)[-groupA]
                }
                twogroup <- c(groupA, groupB)
                lab <- c(rep(1, nA), rep(0, nB))
                t <- mt.teststat(expr.mat[, twogroup], lab, test = method)
                entity <- c(entity1, entity2)
            }
    }
    else {
        t <- NULL
        entity <- NULL
        i <- which(type == "CCR")
        for (j in (1:n.type)[-i]) {
            if (length(which(entity.labels == i)) > minSampleNum) 
                if (length(which(entity.labels == j)) > minSampleNum) 
                  if (type[j] != "Censored") {
                    #cat(".")
                    thisentity1 <- paste(dataset, type[i], sep = "")
                    thisentity2 <- paste(dataset, type[j], sep = "")
                    groupA <- which(entity.labels == i)
                    groupB <- which(entity.labels == j)
                    twogroup <- c(groupA, groupB)
                    nA <- length(groupA)
                    nB <- length(groupB)
                    if (random) {
                      groupA <- sample(twogroup, nA)
                      groupB <- twogroup[!twogroup %in% groupA]
                      twogroup <- c(groupA, groupB)
                    }
                    lab <- c(rep(1, nA), rep(0, nB))
                    thist <- mt.teststat(expr.mat[, twogroup], 
                      lab, test = method)
                    t <- cbind(t, thist)
                    thisentity <- rbind(thisentity1, thisentity2)
                    entity <- cbind(entity, thisentity)
                  }
        }
    }
    if (!random) 
        return(list(entity = entity, stat = t))
    else return(stat = t)
}