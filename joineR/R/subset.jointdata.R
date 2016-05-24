subset.jointdata <- function (x, subj.subset, ...) 
{
    id <- subj.subset
    re <- x
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
    class(re) <- c("jointdata", "list")
    return(re)
}
