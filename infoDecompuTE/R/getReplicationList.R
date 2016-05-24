# get a list of replication for each treatment
getReplicationList <- 
function(design.df, trtCols) {
    
    if (length(trtCols) == 1) {
        
        return(as.matrix(table(design.df[, trtCols])))
    } else if (any(grepl(":", trtCols))) {
        
        level = t(sapply(strsplit(levels(interaction(design.df[, unique(unlist(strsplit(trtCols, "\\:")))])), "\\."), rbind))
        colnames(level) = unique(unlist(strsplit(trtCols, "\\:")))
        
        inter = trtCols[grepl(":", trtCols)]
        
        for (i in 1:length(inter)) {
            level = cbind(level, apply(level[, unique(unlist(strsplit(inter[i], "\\:")))], 1, function(x) paste(x, collapse = ".")))
        }
        
        trtColsList = lapply(strsplit(trtCols, "\\:"), function(x) design.df[, x])
        repList = sapply(trtColsList, function(y) if (is.factor(y)) {
            table(y)
        } else {
            table(apply(y, 1, function(x) paste(x, collapse = ".")))
        })
        
        repList = lapply(repList, function(y) apply(level, 2, function(x) y[x]))
        repList = c(repList, recursive = TRUE)
        repList = matrix(repList[which(!is.na(repList))], nrow = nrow(level))
        
        levelList = sapply(trtColsList, function(y) if (is.factor(y)) {
            nlevels(y)
        } else {
            nlevels(as.factor(apply(y, 1, function(x) paste(x, collapse = "."))))
        })/nlevels(interaction(design.df[, unique(unlist(strsplit(trtCols, "\\:")))]))
        
        repList = repList %*% diag(levelList)
        
        return(repList)
        
    } else {
        level = t(sapply(strsplit(levels(interaction(design.df[, trtCols])), "\\."), rbind))
        colnames(level) = trtCols
        
        repList = apply(design.df[, trtCols], 2, function(x) list(table(x)))
        
        repList = lapply(repList, function(y) apply(level, 2, function(x) y[[1]][x]))
        repList = c(repList, recursive = TRUE)
        repList = matrix(repList[which(!is.na(repList))], nrow = nrow(level))
        
        levelList = apply(design.df[, trtCols], 2, function(x) nlevels(as.factor(x)))/nlevels(interaction(design.df[, trtCols]))
        
        repList = repList %*% diag(levelList)
        
        return(repList)
        
    }
} 
