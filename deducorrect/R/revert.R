
# revert corrections.
# NOTE: it is assumed that actual corrections have taken place. The re
# d     : deducorrect object
# rows  : logical or integer vector indexing records to be reverted
#
revert <- function(d, rows){
    if (missing(rows)) rows <- 1:nrow(d$corrected)
    if (length(rows)==0) return(d)
    if ( is.logical(rows) )  rows <- which(rows)
    status <- d$status
    rows <- rows[status$status[rows] %in% c('corrected','partial')]
    corr <- d$corrections
    cord <- d$corrected
    irws <- corr$row %in% rows
    irow <- corr$row[irws]

    cls <- as.character(corr[irws,'variable'])
    vars <- unique(cls)
    icol <- match(cls,vars)
    A <- as.matrix(cord[vars])

    A[cbind(irow,icol)] <- corr[irws,'old']
    cord[vars] <- A[,vars,drop=FALSE]

    status[rows,'status'] <- "invalid"
    if ( !is.null(status$imputed) ) status[rows,'imputed'] <- 0

    newdeducorrect(
        corrected = cord,
        corrections=corr[!irws,,drop=FALSE],
        status=status,
        Call = d$call
    )
}



correctAndRevert <- function(fun, E, dat, ...){
    v1 <- violatedEdits(E, dat)
    d1 <- fun(E$num, dat, Call=sys.call(-2), ...)
    v2 <- violatedEdits(E, d1$corrected)
    k <- apply(!v1 & v2, 1, any)
    k <- k[!is.na(k)]
    if ( any(k) ) d1 <- revert(d1,rows=k)
    d1
}





