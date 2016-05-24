`sortMatrix` <- 
function(submatrix) {
    
    if (all(is.na(submatrix)) | all(!is.na(submatrix)) | !is.matrix(submatrix)) {
        submatrix
        }
    else {
        submatrix <- submatrix[order(!is.na(submatrix[, ncol(submatrix)])), ]
        lastcol <- ncol(submatrix)
        nas <- is.na(submatrix[, lastcol])
        if (any(nas)) {
            submatrix[which(!nas), -lastcol] <- Recall(submatrix[which(!nas), -lastcol])
            submatrix[which(nas) , -lastcol] <- Recall(submatrix[which(nas) , -lastcol])
            }
        submatrix
        }
    }

