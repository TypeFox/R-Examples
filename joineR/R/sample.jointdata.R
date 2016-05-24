sample.jointdata <-
function (object, size, replace = FALSE) 
{
    origid <- object$subject
    id <- sample(object$subject, size = size, replace = replace)
    re <- object
    if (replace) {
        re$subject <- re$subject[id]
        if (is.data.frame(re$longitudinal)) {
            nn <- diff(match(object$subject, object$longitudinal[[object$subj.col]]))
            nn[length(nn) + 1] <- length(object$longitudinal[, 
                1]) - sum(nn)
		mid <- match(id, origid)
            geti <- match(id, re$longitudinal[[re$subj.col]])
		getirep <- rep(geti, nn[mid])
            ninc <- sequence(nn[mid]) - 1
            lid <- getirep + ninc
            re$longitudinal <- re$longitudinal[lid, ]
            re$longitudinal[[re$subj.col]] <- rep(1:length(id), 
                nn[mid])
            row.names(re$longitudinal) <- 1:(dim(re$longitudinal)[1])
        }
        if (is.data.frame(re$survival)) {
            re$survival <- re$survival[mid, ]
            re$survival[, 1] <- 1:length(id)
            row.names(re$survival) <- 1:(dim(re$survival)[1])
        }
        if (is.data.frame(re$baseline)) {
            re$baseline <- re$baseline[mid, ]
            re$baseline[[re$subj.col]] <- 1:length(id)
            row.names(re$baseline) <- 1:(dim(re$baseline)[1])
        }
    }
    else {
        re$subject <- re$subject[re$subject %in% id]
        if (is.data.frame(re$longitudinal)) {
            re$longitudinal <- re$longitudinal[re$longitudinal[[re$subj.col]] %in% 
                id, ]
            row.names(re$longitudinal) <- 1:(dim(re$longitudinal)[1])
        }
        if (is.data.frame(re$survival)) {
            re$survival <- re$survival[re$survival[[re$subj.col]] %in% 
                id, ]
            row.names(re$survival) <- 1:(dim(re$survival)[1])
        }
        if (is.data.frame(re$baseline)) {
            re$baseline <- re$baseline[re$baseline[[re$subj.col]] %in% 
                id, ]
            row.names(re$baseline) <- 1:(dim(re$baseline)[1])
        }
    }
    class(re) <- c("jointdata", "list")
    return(re)
}
