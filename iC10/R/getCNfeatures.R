getCNfeatures <-
function(CN, Probes, Map, by.feat, ref, Synonyms) {

    if (by.feat=="gene") {
        res <- matrix(NA, nrow=length(Probes), ncol=ncol(CN))
        counter <- 1
        no.na <- 0
        for (i in Probes) {
            ids <- which(tolower(rownames(CN)) == tolower(i))
            if  (length(ids) > 0) {
                vars <- apply(CN[ids,,drop=F], 1, IQR, na.rm=TRUE)
                selected <- which.max(vars)
                res[counter,] <- unlist(CN[ids,,drop=F][selected,])
                no.na <- no.na + 1
            }
            counter <- counter + 1
        }
        if (no.na < length(Probes)) {
            na.ids <- which(apply(res, 1, function(x) mean(is.na(x)))==1)
            for (i in na.ids) {
                one.synon.found <- FALSE
                for (j in Synonyms[which(names(Synonyms)==Probes[i])]) {
                    if (!one.synon.found) {
                    ids <- which(tolower(rownames(CN)) == tolower(j))
                    if  (length(ids) > 0) {
                        vars <- apply(CN[ids,,drop=F], 1, IQR, na.rm=TRUE)
                        selected <- which.max(vars)
                        res[i,] <- unlist(CN[ids,,drop=F][selected,])
                        no.na <- no.na + 1
                    }
                    one.synon.found <- TRUE
                }
                }
            }
        }
        rownames(res) <- Probes
        colnames(res) <- colnames(CN)
        new.D1 <- res
        cat("Found ", no.na, " out of ", length(Probes), " copy number features\n")
    } else {
        Map <- Map[which(Map$Probe_ID %in% Probes),]
        new.D1 <- data.frame(ID=unique(CN$ID))

    for (i in Probes) {
        id <- which(Map$Probe_ID == i)[1] ## Some probes have two annotations
        sub.CN <- CN[which(CN$chromosome_name==Map[id, paste('chromosome_name', ref, sep="_")] &
                           ((CN$loc.start>=Map[id, paste('start_position', ref, sep="_")] &
                           CN$loc.start<=Map[id, paste('end_position', ref, sep="_")]) |
                           (CN$loc.start<=Map[id, paste('start_position', ref, sep="_")] &
                            CN$loc.end>=Map[id, paste('start_position', ref, sep="_")]))),]
        feo <- sapply(split(sub.CN, sub.CN$ID), function(x) {
            if (nrow(x)>0) res <- x$seg.mean[which.max(abs(x$seg.mean))]
            else res <- NA
        })
        feo <- data.frame(ID=names(feo), feo)
        colnames(feo)[2] <- as.character(Map[id,'Probe_ID'])
        new.D1 <- merge(new.D1, feo, all.x=TRUE)
    }
    ## Impute NAs
        if (sum(is.na(new.D1))>0) {
        for (i in 1:nrow(new.D1)) {
        for (j in 1:ncol(new.D1)) {
            if (is.na(new.D1[i,j])) {
                id <- which(Map$Probe_ID == colnames(new.D1)[j])[1] ## Some probes have two annotations
                sub.CN <- CN[which(CN$ID == new.D1[i,1] &
                                   CN$chromosome_name==Map[id, paste('chromosome_name', ref, sep="_")] &
                                   (CN$loc.start>=Map[id, paste('start_position', ref, sep="_")])),][1,]
                new.D1[i,j] <- sub.CN$seg.mean
            }
        }
    }
    }
        rownames(new.D1) <- new.D1$ID
        new.D1$ID <- NULL
        new.D1 <- t(new.D1)
    }
    CN <- new.D1
    CN
}
