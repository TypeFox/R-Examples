getExpfeatures <-
function(Exp, Probes, Synonyms, by.feat) {
    if (by.feat=="probe") {
        Exp <- Exp[which(rownames(Exp) %in% Probes),]
        Exp <- Exp[match(Probes, rownames(Exp)),]
    } else {
        DS.gene <- matrix(NA, ncol=ncol(Exp), nrow=length(Probes))
        counter <- 1
        no.na <- 0
        for (i in Probes) {
            ids <- which(tolower(rownames(Exp)) == tolower(i))
            if  (length(ids) > 0) {
                vars <- apply(Exp[ids,,drop=F], 1, IQR, na.rm=TRUE)
                if (mean(1 * (is.na(vars)))<1) {
                    selected <- which.max(vars)
                    DS.gene[counter,] <- unlist(Exp[ids,,drop=F][selected,])
                    no.na <- no.na + 1
                }
            }
            counter <- counter + 1
        }
        if (no.na < length(Probes)) {
            na.ids <- which(apply(DS.gene, 1, function(x) mean(is.na(x)))==1)
            for (i in na.ids) {
                one.synon.found <- FALSE
                for (j in Synonyms[which(names(Synonyms)==Probes[i])]) {
                    if (!one.synon.found) {
                        ids <- which(tolower(rownames(Exp)) == tolower(j))
                        if  (length(ids) > 0) {
                            vars <- apply(Exp[ids,,drop=F], 1, IQR, na.rm=TRUE)
                            if (mean(1 * (is.na(vars)))<1) {
                                selected <- which.max(vars)
                                DS.gene[i,] <- unlist(Exp[ids,,drop=F][selected,])
                                no.na <- no.na + 1
                            }
                            one.synon.found <- TRUE
                        }
                    }
                }
            }
        }
        rownames(DS.gene) <- Probes
        colnames(DS.gene) <- colnames(Exp)
        Exp <- DS.gene
        cat("Found ", no.na, " out of ", length(Probes), " Exp features\n")
    }
    as.matrix(Exp)
}
