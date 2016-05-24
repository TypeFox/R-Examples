`rowDominance` <-
function(mtrx) {
    
    rownames(mtrx) <- seq(nrow(mtrx))
    mtrx.copy <- mtrx
    sums <- rowSums(mtrx)
    mtrx <- mtrx[order(sums, decreasing=TRUE), , drop=FALSE]
    
    sums <- sort(sums, decreasing=TRUE)
    line.no <- 1
    while(line.no < nrow(mtrx)) {
        less <- sums < sums[line.no]
        if (any(less)) {
            aa <- apply(mtrx[less, , drop=FALSE], 1, function(x) {all(mtrx[line.no, x])})
            mtrx <- rbind(mtrx[!less, , drop=FALSE], mtrx[less, , drop=FALSE][!aa, , drop=FALSE])
            sums <- rowSums(mtrx)
            line.no <- line.no + 1
        }
        else {
            break
        }
    }
    
    return(match(rownames(mtrx), rownames(mtrx.copy)))
}

