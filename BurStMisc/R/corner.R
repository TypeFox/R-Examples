"corner" <-
function (x, corner = "tlffff", n = rep(6, length(corner))) 
{
    if (length(corner) != 1) 
        stop("corner must be a single string")
    dx <- dim(x)
    ldx <- length(dx)
    if (ldx < 2) {
        if (substring(corner, 1, 1) == "t") 
            return(head(x, n = n[1]))
        else return(tail(x, n = n[1]))
    }
    clen <- nchar(corner)
    if (clen < ldx) {
        if (clen == 0) {
            corner <- paste("tl", paste(rep("f", ldx - 2), collapse = ""), 
                sep = "")
        }
        else if (clen == 1) 
            corner <- paste(corner, "l", paste(rep("f", ldx - 
                1), collapse = ""), sep = "")
        else corner <- paste(corner, paste(rep("f", ldx - clen), 
            collapse = ""), sep = "")
    }
    corner <- substring(corner, 1:ldx, 1:ldx)
    n <- rep(n, length = ldx)
    if (corner[1] == "t") 
        rsub <- seq(length = min(n[1], dx[1]))
    else rsub <- seq(to = dx[1], length = min(n[1], dx[1]))
    if (corner[2] == "l") 
        csub <- seq(length = min(n[2], dx[2]))
    else csub <- seq(to = dx[2], length = min(n[2], dx[2]))
    if (ldx == 2) {
        return(x[rsub, csub, drop = FALSE])
    }
    subv <- vector("list", ldx + 1)
    subv[[1]] <- rsub
    subv[[2]] <- csub
    for (i in 3:ldx) {
        if (corner[i] == "f") 
            subv[[i]] <- seq(length = min(n[i], dx[i]))
        else subv[[i]] <- seq(to = dx[i], length = min(n[i], 
            dx[i]))
    }
    names(subv) <- c(rep("", ldx), "drop")
    subv[[ldx + 1]] <- FALSE
    do.call("[", c(list(x), subv))
}
