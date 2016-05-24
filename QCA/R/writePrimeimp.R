 ## when the matrix cannot be further minimized, this function writes the prime implicants
 ## as the name of the conditions (columns), collapsed together in a single string

`writePrimeimp` <- 
function(idx, collapse="*", uplow=FALSE, use.tilde=FALSE) {
    
    if (use.tilde) {
        uplow <- FALSE
    }
    
    idx <- as.data.frame(idx)
    
    for (i in seq(ncol(idx))) {
        if (uplow) {
            conditions <- c(tolower(colnames(idx)[i]), toupper(colnames(idx)[i]))
        }
        else if (use.tilde) {
            conditions <- c(paste("~", toupper(colnames(idx)[i]), sep=""), toupper(colnames(idx)[i]))
        }
        else {
            conditions <- paste(colnames(idx)[i], "{", seq(max(idx[, i])) - 1, "}", sep="")
        }
        idx[idx[, i] != 0, i] <- conditions[idx[idx[, i] != 0, i]]
    }
    
    return(apply(idx, 1, function(x) {
        paste(x[x != 0], collapse=collapse)
    }))
}

