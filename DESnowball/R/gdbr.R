gdbr <- function(y,D) {
    A = (-1/2) * D^2
    n <- length(y)
    I = diag(1, n)
    A.cen = scale(A, scale=FALSE)
    G = A.cen - rowMeans(A.cen)

    gdbr.stat = gdbr_fstat(y, G)
    gdbr.stat
}

gdbr_fstat <- function(casecon, G) {
    y.new = casecon - mean(casecon)
    I = diag(1, length(casecon))
    H = y.new %*% solve((t(y.new) %*% y.new)) %*% t(y.new)
    Fstat.num = sum(diag(H %*% G %*% H))
    Fstat.denom = sum(diag((I - H) %*% G %*% t(I - H)))
    Fstat = Fstat.num / Fstat.denom	
    Fstat
}


