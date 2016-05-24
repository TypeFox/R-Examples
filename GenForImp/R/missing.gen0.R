missing.gen0 <-
function (mat, nummiss)
{
    n <- nrow(mat)
    p <- ncol(mat)
    mmiss <- mat
    posmiss <- sample(n * p, nummiss, replace = FALSE)-1
    # row index
    rowpos <- posmiss%%n + 1
    # column index
    colpos <- posmiss%/%n+1
    # (row index, column index)
    ind <- array(c(rowpos,colpos),dim=c(nummiss,2))
    indmat <- matrix(FALSE,n,p)
    indmat[ind] <- TRUE
    mmiss[indmat] <- NA
    mmiss
}
