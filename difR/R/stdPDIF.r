# DIF: STANDARDIZATION
stdPDIF<-function (data, member, anchor = 1:ncol(data), stdWeight="focal") 
{
    resPDIF <- resAlpha<-NULL
    for (item in 1:ncol(data)) {
        data2 <- data[, anchor]
        if (sum(anchor == item) == 0) 
            data2 <- cbind(data2, data[, item])
        xj <- rowSums(data2, na.rm=TRUE)
        scores <- sort(unique(xj))
        ind <- 1:nrow(data)
        prov <- NULL
        for (j in 1:length(scores)) {
            Prs <- length(ind[xj == scores[j] & member == 0 & 
                data[, item] == 1])/length(ind[xj == scores[j] & 
                member == 0])
            Pfs <- length(ind[xj == scores[j] & member == 1 & 
                data[, item] == 1])/length(ind[xj == scores[j] & 
                member == 1])
            Ks <- switch(stdWeight,
			focal=length(ind[xj == scores[j] & member == 1]),
			reference=length(ind[xj == scores[j] & member == 0]),
			total=length(ind[xj == scores[j]]))
		Nfs <- length(ind[xj == scores[j] & member == 1])
            prov <- rbind(prov, c(scores[j], Prs, Pfs, Ks, Nfs))
        }
        stdNum <- stdDen <- pfNum <- prNum <- pDen <- 0
        for (i in 1:nrow(prov)) {
            if (is.na(prov[i, 2]) == FALSE & is.na(prov[i, 3]) == 
                FALSE) {
                stdNum <- stdNum + prov[i, 4] * (prov[i, 3] - prov[i, 
                  2])
                stdDen <- stdDen + prov[i, 4]
		    pfNum<-pfNum+ prov[i,5]*prov[i,3]
		    pDen<-pDen+ prov[i,5]
		    prNum<-prNum+ prov[i,5]*prov[i,2]
            }
        }
	  Pf<-pfNum/pDen
	  Pr<-prNum/pDen
        resPDIF[item] <- stdNum/stdDen
	  resAlpha[item] <- (Pr/(1-Pr))/(Pf/(1-Pf)) 
    }
    return(list(resStd=resPDIF,resAlpha=resAlpha))
}
