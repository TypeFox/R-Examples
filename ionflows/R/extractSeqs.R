extractSeqs <-
function(bed.table, genom) {
    n <- nrow(bed.table)
    seqs <- rep(NA, n)
    cat("Extracting sequences (. = 10 amplicons):\n")
    for(i in 1:n) {
        seqs[i] <- as.character(getSeq(genom, bed.table[i, "chrom"], bed.table[i, "chromStart"] + 1, bed.table[i, "chromEnd"]))
        if(i%%10 == 0) cat(".")
        if(i%%500 == 0 || i == n) cat(i, "/" ,n ,"\n", sep="")
    }
    return(seqs)
}
