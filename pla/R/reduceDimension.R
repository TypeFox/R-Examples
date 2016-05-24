.reduceDimensionSub <- function(A, i, j = i+1) {
    dms <- dim(A)
    nms <- dimnames(A)
    nms[[i]] <- paste(rep(nms[[j]], rep(dms[i], dms[j])),
                      rep(nms[[i]], dms[j]), sep = ":")
    dms[i] <- dms[i] * dms[j]
    dim(A) <- dms[-(j)]
    dimnames(A) <- nms[-(j)]
    return(A)
}

.reduceDimension <- function(Array, namesPermuted = dimnames(Array),
                             selectFun = function (array) NULL) {
    if (length(dim(Array)) > 2)
        Assay <- .reduceDimensionSub(Array, length(dim(Array))-1)
    else
        Assay <- Array
    ldim <- length(dim(Assay))
    nmsD <- namesPermuted[ldim][[1]]
    nmsS <- namesPermuted[ldim+1][[1]]
    if (!is.null(body(selectFun)))
        if ((dim(Assay)[ldim] < length(nmsD)) |
            (dim(Assay)[ldim] < length(nmsS)))
            warning(paste("Unable to correct 'dose' and 'sample' labels\n",
                          "of levels when doing selection on these factors"),
                    call. = FALSE)
    if (dim(Assay)[ldim] == length(nmsD)) {
        if (dim(Assay)[ldim] == length(nmsS)) {
            idxS <- unlist(lapply(diff(c(match(unique(nmsS), nmsS),
                                         length(nmsS)+1)),
                                  FUN = function(x) 1:x))
            nmsSD <- paste(nmsS, idxS, nmsD, sep = ":")
            dimnames(Assay)[[ldim]] <- nmsSD
        } else {
            dimnames(Assay)[[ldim]] <- nmsD
        }
    }
    return(Assay)
}
