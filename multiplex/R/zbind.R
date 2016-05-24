zbind <-
function (...) 
{
    argl <- list(...)
    ifelse(isTRUE(dim(argl[[1]])[3] > 1) == TRUE, pvt <- argl[[1]][, 
        , 1], pvt <- argl[[1]])
    tmp <- data.frame(matrix(ncol = (dim(pvt)[1] * dim(pvt)[2]), 
        nrow = 0))
    for (i in 1:length(argl)) {
        if (isTRUE(is.array(argl[[i]]) == TRUE) == FALSE) 
            stop("(part of) 'x' is not an array or matrix")
        if (isTRUE(dim(argl[[i]])[3] > 1) == TRUE) {
            for (j in 1:dim(argl[[i]])[3]) {
                if (isTRUE(dim(argl[[i]][, , j])[1] != dim(pvt)[1]) == 
                  TRUE) 
                  stop("Dimensions 'x', 'y' must be equal")
                tmp[(nrow(tmp) + 1), ] <- as.vector(argl[[i]][, 
                  , j])
            }
            rm(j)
        }
        else if (isTRUE(is.na(dim(argl[[i]])[3]) == TRUE) == 
            TRUE | isTRUE(dim(argl[[i]])[3] == 1) == TRUE) {
            if (isTRUE(dim(argl[[i]])[1] != dim(pvt)[1]) == TRUE) 
                stop("Dimensions 'x', 'y' must be equal")
            tmp[(nrow(tmp) + 1), ] <- as.vector(argl[[i]])
        }
    }
    rm(i)
    arr <- array(dim = c(dim(pvt)[1], dim(pvt)[2], nrow(tmp)))
    for (i in 1:nrow(tmp)) {
        arr[, , i][1:(dim(pvt)[1] * dim(pvt)[2])] <- as.numeric(tmp[i, 
            ])
    }
    rm(i)
    if (isTRUE(is.null(dimnames(arr)) == TRUE) == TRUE) {
        if (isTRUE(is.null(dimnames(pvt)) == TRUE) == FALSE) 
            dimnames(arr)[[1]] <- dimnames(arr)[[2]] <- dimnames(pvt)[[1]]
    }
    lbs <- vector()
    for (i in 1:length(argl)) {
        if (isTRUE(dim(argl[[i]])[3] > 1) == TRUE) {
            ifelse(is.null(dimnames(argl[[i]])[[3]]) == FALSE, 
                lbs <- append(lbs, dimnames(argl[[i]])[[3]]), 
                lbs <- append(lbs, (length(lbs) + 1:dim(argl[[i]])[3])))
        }
        else if (isTRUE(dim(argl[[i]])[3] > 1) == FALSE) {
            lbs <- append(lbs, (length(lbs) + 1))
        }
    }
    rm(i)
    if (isTRUE(length(lbs) == dim(arr)[3]) == TRUE) 
        dimnames(arr)[[3]] <- as.list(lbs)
    arr
}
