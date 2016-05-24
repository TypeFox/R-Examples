Isoqval <- 
function (delta, allfdr, qqstat, stat) 
{
    qstat <- switch(stat, E2 = qqstat[[1]], Williams = qqstat[[3]], 
        Marcus = qqstat[[5]], M = qqstat[[7]], ModifM = qqstat[[9]])
    obs.stat <- qstat[order(qstat[, 4]), 1]
    if (sum(qstat[qstat[, 1] < 0, 3] < -delta) > 0) {
        low <- max(qstat[qstat[, 1] < 0, 1][qstat[qstat[, 1] < 
            0, 3] < -delta])
    }
    else {
        low <- -Inf
    }
    if (low != "-Inf") {
        which.low <- which(qstat[qstat[, 1] < 0, 1] == low)
        low.point <- which.low[length(which.low)]
    }
    if (low == "-Inf") {
        low <- min(qstat[, 1])
        low.point <- 0
    }
    if (sum(qstat[qstat[, 1] > 0, 3] > delta) > 0) {
        up <- min(qstat[qstat[, 1] > 0, 1][qstat[qstat[, 1] > 
            0, 3] > delta])
    }
    else {
        up <- Inf
    }
    if (up == "Inf") {
        up <- max(qstat[qstat[, 1] > 0, 1])
        up.point <- nrow(qstat) + 1
    }
    if (up != "Inf") {
        which.up <- which(qstat[qstat[, 1] > 0, 1] == up)
        up.point <- sum(qstat[, 1] <= 0) + which.up[1]
    }
    dtable <- allfdr
    k <- nrow(dtable) - 1
    nnew <- matrix(0, nrow(qstat), 2)
    test <- abs(qstat[, 3]) < dtable[1, 1]
    nnew[test, 1] <- qstat[test, 4]
    nnew[test, 2] <- dtable[1, 5]
    test2 <- abs(qstat[, 3]) >= dtable[k + 1, 1]
    nnew[test2, 1] <- qstat[test2, 4]
    nnew[test2, 2] <- dtable[k + 1, 5]
    for (i in k:1) {
        j <- i + 1
        test1 <- dtable[i, 1] <= abs(qstat[, 3]) & abs(qstat[, 
            3]) < dtable[j, 1]
        nnew[test1, 1] <- qstat[test1, 4]
        nnew[test1, 2] <- dtable[i, 5]
    }
    nnnew <- cbind(nnew, qstat[, 3], nnew[, 2])
    fdr <- NULL
    for (i in k:1) {
        j <- i + 1
        test2 <- dtable[i, 1] <= delta & delta < dtable[j, 1]
        if (sum(test2) > 0) {
            fdr <- dtable[i, 5]
        }
        test3 <- dtable[i, 1] < delta & delta <= dtable[j, 1]
        if (sum(test3) > 0) {
            fdr <- dtable[j, 5]
        }
    }
    n1 <- n2 <- NULL
    if (low.point == 0 & is.numeric(fdr)) {
        n1 <- NULL
    }
    if (low.point == 1 & is.numeric(fdr)) {
        n1 <- matrix(nnnew[1, ], 1, 4)
    }
    if (low.point > 1 & is.numeric(fdr)) {
        n1 <- nnnew[1:low.point, ]
        n1[n1[, 3] > -delta, 4] <- fdr
    }
    if (up.point == nrow(qstat) + 1 & is.numeric(fdr)) {
        n2 <- NULL
    }
    if (up.point == nrow(qstat) & is.numeric(fdr)) {
        n2 <- matrix(nnnew[nrow(qstat), ], 1, 4)
    }
    if (up.point < nrow(qstat) & is.numeric(fdr)) {
        n2 <- nnnew[up.point:nrow(qstat), ]
        n2[n2[, 3] < delta, 4] <- fdr
    }
    res <- sign.list <- res1 <- sign.list1 <- NULL
    m1 <- low.point + 1
    m2 <- up.point - 1

    cols <- c(1, 4)
    if (is.numeric(fdr) & is.numeric(n1)) {
        n1[n1[, 4] > fdr, 4] <- fdr
    }
    if (is.numeric(fdr) & is.numeric(n2)) {
        n2[n2[, 4] > fdr, 4] <- fdr
    }
   
       

 if (is.numeric(fdr) & is.numeric(n1) & is.numeric(n2)) {
         if (m2 >= m1) { 
       res <- rbind(n1[, cols], nnnew[m1:m2, cols], n2[, 
                cols]) } else
        res <- rbind(n1[, cols], n2[, cols])
            
        res1 <- cbind(res[, 1], obs.stat[res[, 1]], res[, 2])
            sign.list <- rbind(n1[, cols], n2[, cols])
            sign.list1 <- cbind(sign.list[, 1], obs.stat[sign.list[, 
                1]], sign.list[, 2])
        }


        if (is.numeric(fdr) & is.numeric(n1) & !is.numeric(n2)) {
           if (m2 >= m1) { res <- rbind(n1[, cols], nnnew[m1:m2, cols])} else
           res <- n1[, cols]
            res1 <- cbind(res[, 1], obs.stat[res[, 1]], res[, 
                2])
            sign.list <- rbind(n1[, cols])
            sign.list1 <- cbind(sign.list[, 1], obs.stat[sign.list[, 
                1]], sign.list[, 2])
        }



        if (is.numeric(fdr) & !is.numeric(n1) & is.numeric(n2)) {
           if (m2 >= m1) {  res <- rbind(nnnew[m1:m2, cols], n2[, cols])} else
               res <- n2[, cols]
            res1 <- cbind(res[, 1], obs.stat[res[, 1]], res[, 
                2])
            sign.list <- rbind(n2[, cols])
            sign.list1 <- cbind(sign.list[, 1], obs.stat[sign.list[, 
                1]], sign.list[, 2])
        }
    
    colnames(res1) <- colnames(sign.list1) <- c("Row.names", 
        "t.stat", "q.val")
    return(list(res1, sign.list1))
}
