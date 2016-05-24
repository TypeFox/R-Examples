zbnd <-
function (a, b) 
{
    if (isTRUE(is.na(dim(a)[3]) == FALSE)) {
        if (isTRUE(dim(b)[3] > 1) == TRUE) {
            nx <- array(a, dim = c(dim(a)[1], dim(a)[2], (dim(a)[3] + 
                dim(b)[3])))
            nx[, , (dim(a)[3] + 1):(dim(nx)[3])] <- b
        }
        else {
            nx <- array(a, dim = c(dim(a)[1], dim(a)[2], (dim(a)[3] + 
                1)))
            nx[, , (dim(a)[3] + 1)] <- b
        }
    }
    else {
        if (isTRUE(dim(b)[3] > 1) == TRUE) {
            nx <- array(a, dim = c(dim(a)[1], dim(a)[2], (1 + 
                dim(b)[3])))
            nx[, , 2:(dim(nx)[3])] <- b
        }
        else {
            nx <- array(a, dim = c(dim(a)[1], dim(a)[2], 2))
            nx[, , 2] <- b
        }
    }
    return(nx)
}
