signed <-
function (P, N = NULL, labels = NULL) 
{
    if (is.array(P) == FALSE) 
        stop("Data must be an array")
    if (is.null(N) == FALSE) {
        if (is.na(dim(N)[3]) == FALSE) {
            N <- N[, , 1]
            warning("Take the 1st dim. in 'N' only.")
        }
    }
    else if (is.null(N) == TRUE) {
        if (is.na(dim(P)[3]) == TRUE) {
            N <- dichot(P, c = max(P) + 1L)
            warning("No negative ties are provided.")
            Po <- P
        }
        else {
            N <- P[, , 2]
            Po <- P <- P[, , 1]
        }
    }
    P <- dichot(P, c = 1L)
    N <- dichot(N, c = 1L)
    sm <- P + N
    if (isTRUE(any((sm) > 1L)) == TRUE) {
        ambs <- which(sm == 2L)
        sm[which(P == 1L)] <- "p"
        sm[which(N == 1L)] <- "n"
        sm[ambs] <- "a"
        sm[which(suppressWarnings(as.numeric(sm) == 0))] <- "o"
    }
    else {
        sm[which(N == 1L)] <- -1L
    }
    if (is.null(labels) == FALSE) {
        ifelse(isTRUE(length(labels) == dim(sm)[1]) == TRUE, 
            NA, labels <- 1:dim(sm)[1])
        rownames(sm) <- colnames(sm) <- labels
    }
    else if (is.null(dimnames(P)[1]) == FALSE) {
        rownames(sm) <- colnames(sm) <- dimnames(P)[[1]]
    }
    else {
        rownames(sm) <- colnames(sm) <- 1:dim(sm)[1]
    }
    val <- levels(factor(sm))
    lst <- list(val = noquote(val[length(val):1]), s = noquote(sm))
    class(lst) <- "Signed"
    return(lst)
}
