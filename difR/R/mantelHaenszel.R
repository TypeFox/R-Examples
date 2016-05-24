mantelHaenszel<-function (data, member, correct = TRUE, exact=FALSE,anchor = 1:ncol(data)) 
{
    res <- resAlpha <- varLambda <- RES<-NULL
    for (item in 1:ncol(data)) {
        data2 <- data[, anchor]
        if (sum(anchor == item) == 0) 
            data2 <- cbind(data2, data[, item])
if (!is.matrix(data2)) data2<-cbind(data2)
        xj <- rowSums(data2, na.rm = TRUE)
        scores <- sort(unique(xj))
        prov <- NULL
        ind <- 1:nrow(data)
        for (j in 1:length(scores)) {
            Aj <- length(ind[xj == scores[j] & member == 0 & 
                data[, item] == 1])
            Bj <- length(ind[xj == scores[j] & member == 0 & 
                data[, item] == 0])
            Cj <- length(ind[xj == scores[j] & member == 1 & 
                data[, item] == 1])
            Dj <- length(ind[xj == scores[j] & member == 1 & 
                data[, item] == 0])
            nrj <- length(ind[xj == scores[j] & member == 0])
            nfj <- length(ind[xj == scores[j] & member == 1])
            m1j <- length(ind[xj == scores[j] & data[, item] == 
                1])
            m0j <- length(ind[xj == scores[j] & data[, item] == 
                0])
            Tj <- length(ind[xj == scores[j]])
if (exact){
if (Tj > 1) prov <- c(prov, c(Aj,Bj, Cj, Dj))
}
else{
            if (Tj > 1) 
                prov <- rbind(prov, c(Aj, nrj * m1j/Tj, (((nrj * 
                  nfj)/Tj) * (m1j/Tj) * (m0j/(Tj - 1))), scores[j], 
                  Bj, Cj, Dj, Tj))
}
        }
if (exact){
tab<-array(prov,c(2,2,length(prov)/4))
pr<-mantelhaen.test(tab,exact=TRUE)
RES<-rbind(RES,c(item,pr$statistic,pr$p.value))

}
else{
        if (correct == TRUE) 
            res[item] <- (abs(sum(prov[, 1] - prov[, 2])) - 0.5)^2/sum(prov[, 
                3])
        else res[item] <- (abs(sum(prov[, 1] - prov[, 2])))^2/sum(prov[, 
            3])
        resAlpha[item] <- sum(prov[, 1] * prov[, 7]/prov[, 8])/sum(prov[, 
            5] * prov[, 6]/prov[, 8])
        varLambda[item] <- sum((prov[, 1] * prov[, 7] + resAlpha[item] * 
            prov[, 5] * prov[, 6]) * (prov[, 1] + prov[, 7] + 
            resAlpha[item] * (prov[, 5] + prov[, 6]))/prov[, 
            8]^2)/(2 * (sum(prov[, 1] * prov[, 7]/prov[, 8]))^2)
    }
}
if (exact) return(list(resMH=RES[,2],Pval=RES[,3]))
else    return(list(resMH = res, resAlpha = resAlpha, varLambda = varLambda))
}
