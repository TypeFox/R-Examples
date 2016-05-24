make.contrasts <-
function (data = data[, loci], allele.chars = letters)
{
    subsets <- function(n, r, v = 1:n) if (r <= 0)
        vector(mode(v), 0)
    else if (r >= n)
        v[1:n]
    else {
        rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), Recall(n -
                                                               1, r, v[-1]))
    }
    data <- data.frame(data)
    loci <- names(data)
    nobs <- dim(data)[1]
    oset <- rep(1,nobs)
    nloci <- length(loci)
    contr.df <- cbind(data, matrix(0, nrow = nobs, ncol = 2 * nloci^2))
    aterms <- paste(allele.chars[1:nloci], sep = "")
    aaterms <- paste(aterms, allele.chars[1:nloci], sep = "")
    k <- nloci + 1
    n1 <- length(aterms) + length(aaterms)
    list.columns <- list(aterms = aterms, aaterms = aaterms)
    names(contr.df)[k:(k + n1 - 1)] <- c(aterms, aaterms)
    if (nloci > 1) {
        all2 <- apply(matrix(subsets(nloci, 2), ncol = 2), 1,
                      function(x) allele.chars[x])
        ab <- apply(all2, 2, function(x) paste(x, collapse = ""))
        sabterms <- paste("s", ab, sep = "")
        qabterms <- paste("q", ab, sep = "")
        abbaabterms <- c(paste(ab, all2[2, ], sep = ""), paste(all2[1,
                                                                    ], ab, sep = ""))
        list.columns <- c(list.columns, list(sabterms = sabterms,
                                             qabterms = qabterms, abbaabterms = abbaabterms))
        n2 <- 4 * length(sabterms)
        names(contr.df)[(k + n1):(k + n1 + n2 - 1)] <- c(sabterms,
                                                         qabterms, abbaabterms)
    }
    for (i in 1:nloci) {
        setup <- decode.genotypes(contr.df[, loci[i]])
        char2 <- allele.chars[i]
        n.mb <- paste(char2, sep = "")
        n.mbb <- paste(char2, char2, sep = "")
        contr.df[, n.mb] <- mb <- setup$ma
        contr.df[, n.mbb] <- mbb <- setup$maa
        oset <- setup$oset * oset
        if (i > 1)
            for (j in 1:(i - 1)) {
                char1 <- allele.chars[j]
                n.ma <- char1
                n.maa <- paste(char1, char1, sep = "")
                ma <- contr.df[, n.ma]
                maa <- contr.df[, n.maa]
                ab <- paste(char1, char2, sep = "")
                nsab <- paste("s", ab, sep = "")
                nqab <- paste("q", ab, sep = "")
                nmabb <- paste(ab, char2, sep = "")
                nmaab <- paste(char1, ab, sep = "")
                contr.df[, nsab] <- as.numeric(ma == 1 & mb ==
                                               1)
                contr.df[, nqab] <- (ma > 0 & mb > 0) + (maa *
                                                         mbb) - ((ma == 1) * (mb == 1))
                contr.df[, nmabb] <- ma * (mb == 2)
                contr.df[, nmaab] <- mb * (ma == 2)
            }
    }
    list(contrasts.df = contr.df, oset=oset, list.columns = list.columns)
}
