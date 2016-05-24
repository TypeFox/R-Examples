isolateLoc <-
function (reg, loc, remove0 = FALSE) 
{
    newE <- reg$E
    newE[is.na(newE)] <- 0
    arr <- matrix(unlist(strsplit(names(newE), "")), byrow = TRUE, 
        nrow = length(names(newE)))
    if (length(loc) < reg$nloc) {
        newE[!apply(as.matrix(arr[, -loc]) != ".", 1, sum) == 
            0] <- 0
    }
    if (remove0) {
        newE <- newE[newE != 0]
    }
    return(newE)
}
