sbs1 <-
function (data, indic, correct0, kvec, method = c("lda", "knn", 
    "rpart")) 
{
    n = dim(data)[1]
    p = dim(data)[2]
    output <- indic
    varia <- 1:(p - 1)
    varia <- varia[indic > 0]
    correct <- rep(0, p - 1)
    mm = 0
    for (m in 1:(p - 1)) {
        if (indic[m] == 1) {
            mm = mm + 1
            temp <- varia
            which <- temp[-mm]
            if (method == "lda") 
                correct[m] <- cv10lda2(data[, c(which, p)])
            else if (method == "knn") 
                correct[m] <- cv10knn2(data[, c(which, p)], 
                  kvec)
            else correct[m] <- cv10rpart2(data[, c(which, p)])
        }
    }
    prov <- correct + runif(p - 1)
    where <- sum((1:(p - 1)) * as.numeric(max(prov) == prov))
    output <- correct[where]/n
    if (output >= correct0) {
        indic[where] <- 1
    }
    else {
        output <- correct0
        where <- NULL
        which1 <- NULL
    }
    which <- rev(which)
    which1 <- where
    indic[where] <- 0
    list(variaelim = which1, indic = indic, correcto = output)
}
