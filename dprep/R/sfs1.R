sfs1 <-
function (data, indic, correcto, kvec, method = c("lda", "knn", 
    "rpart")) 
{
    n = dim(data)[1]
    p = dim(data)[2]
    output <- indic
    varia <- 1:(p - 1)
    varia <- varia[indic > 0]
    correct <- rep(0, p - 1)
    for (m in 1:(p - 1)) {
        if (indic[m] == 0) {
            which <- c(m, varia, p)
            if (method == "lda") 
                correct[m] <- cv10lda2(data[, which])
            else if (method == "knn") 
                correct[m] <- cv10knn2(data[, which], kvec)
            else correct[m] <- cv10rpart2(data[, which])
        }
    }
    prov <- correct + runif(p - 1)
    where <- which(max(prov) == prov)
    output <- correct[where]/n
    if (output > correcto) {
        indic[where] <- 1
    }
    list(indic = indic, varselec = where, accuracy = output)
}
