write.srt <-
function (x, file = NULL, sep = "\t", header = TRUE) 
{
    if (isTRUE(is.na(dim(x)[3]) == TRUE) == FALSE) {
        if (isTRUE(is.null(dimnames(x)[[3]])) == TRUE) 
            dimnames(x)[[3]] <- 1:dim(x)[3]
    }
    if (isTRUE(is.null(dimnames(x)[[1]])) == TRUE) 
        dimnames(x)[[1]] <- dimnames(x)[[2]] <- 1:dim(x)[1]
    if (isTRUE(dim(x)[3] == 1) == TRUE) 
        x <- x[, , 1]
    if (isTRUE(is.na(dim(x)[3]) == TRUE) == TRUE) {
        if (header) {
            cat(paste("Sender", "Receiver", "Ties", sep = sep), 
                file = file, sep = "\n", append = TRUE)
        }
        tmp <- transf(x, lb2lb = TRUE, labels = dimnames(x)[[1]])
        for (i in 1:length(tmp)) {
            cat(paste(strsplit(tmp[i], ", ")[[1]][1], strsplit(tmp[i], 
                ", ")[[1]][2], "1", sep = sep), file = file, 
                sep = "\n", append = TRUE)
        }
        rm(i)
    }
    else if (isTRUE(is.na(dim(x)[3]) == TRUE) == FALSE) {
        ts <- rep(0L, dim(x)[3])
        DF <- data.frame(matrix(ncol = (dim(x)[3] + 2L), nrow = 0L))
        for (k in 1:dim(x)[3]) {
            tmp <- transf(x[, , k], lb2lb = TRUE, labels = dimnames(x)[[1]])
            for (i in 1:length(tmp)) {
                ts[k] <- 1L
                DF[(nrow(DF) + 1L), ] <- c(strsplit(tmp[i], ", ")[[1]][1], 
                  strsplit(tmp[i], ", ")[[1]][2], ts)
                ts <- rep(0L, dim(x)[3])
            }
            rm(i)
        }
        rm(k)
        isol <- dimnames(x)[[1]][which(!(dimnames(x)[[1]] %in% 
            dimnames(rm.isol(x))[[1]]))]
        if (isTRUE(length(isol) > 0L) == TRUE) {
            for (i in 1:length(isol)) {
                DF[(nrow(DF) + 1L), ] <- c(isol[i], isol[i], 
                  ts)
            }
            rm(i)
        }
        colnames(DF) <- c("Sender", "Receiver", dimnames(x)[[3]])
        if (header) {
            utils::write.table(DF, file = file, quote = FALSE, 
                sep = sep, row.names = FALSE, col.names = TRUE)
        }
        else if (!header) {
            utils::write.table(DF, file = file, quote = FALSE, 
                sep = sep, row.names = FALSE, col.names = FALSE)
        }
    }
}
