fltr <-
function (x, PO, rclos = TRUE, ideal = FALSE) 
{
    if (is.null(dimnames(PO)[[1]]) == TRUE) 
        stop("Dimnames in 'PO' are NULL")
    if (isTRUE(is.character(x) == TRUE) == TRUE) {
        lbs <- dimnames(PO)[[1]]
        tmp <- jnt(unlist(strsplit(lbs, "} {", fixed = TRUE)), 
            prsep = ", ")
        tmp <- sub("{", "", dhc(tmp, prsep = ", "), fixed = TRUE)
        tmp <- sub("}", "", dhc(tmp, prsep = ", "), fixed = TRUE)
        if (isTRUE(length(tmp) != length(unique(tmp))) == TRUE) 
            stop("'PO' must be in a reduced form.")
        flg <- FALSE
        for (i in 1:length(lbs)) {
            tmp <- jnt(unlist(strsplit(lbs[i], "} {", fixed = TRUE)), 
                prsep = ", ")
            tmp <- sub("{", "", dhc(tmp, prsep = ", "), fixed = TRUE)
            tmp <- sub("}", "", dhc(tmp, prsep = ", "), fixed = TRUE)
            if (isTRUE(x %in% tmp) == TRUE) {
                X <- i
                flg <- TRUE
                break
            }
            else {
                NA
            }
        }
        rm(i)
        if (isTRUE(flg == FALSE) == TRUE) 
            stop("'x' is not part of the given partial order.")
    }
    else {
        if (isTRUE(x > nrow(PO) | x <= 0L) == TRUE) 
            stop("'x' is either non-positive or greater than the size of the given partial order.")
        X <- as.integer(x)
    }
    ifelse(isTRUE(ideal == TRUE) == TRUE, po <- t(PO), po <- PO)
    pfl <- vector()
    for (i in 1:nrow(po)) {
        if (isTRUE(po[X, i] == 1L) == TRUE && isTRUE(po[i, X] == 
            0L) == TRUE) {
            pfl <- append(pfl, i)
        }
        else {
            NA
        }
    }
    rm(i)
    if (rclos) {
        pfl <- append(X, pfl)
    }
    if (isTRUE(length(pfl) > 0L) == TRUE) {
        pfll <- as.list(dimnames(PO)[[1]][pfl])
        attr(pfll, "names") <- pfl
    }
    else {
        pfll <- NULL
    }
    return(pfll)
}
