dematrify <- function (taxa, filename, sep = ",", thresh = 0) 
{
    tmp <- which(taxa > thresh, arr.ind = TRUE)
    samples <- row.names(tmp)
    species <- names(taxa)[tmp[, 2]]
    abund <- taxa[tmp]
    ord <- order(tmp[, 1], tmp[, 2])
    result <- data.frame(samples[ord], species[ord], abund[ord])
    names(result) <- c("sample", "species", "abundance")
    attr(result,'call') <- match.call()
    attr(result,'taxa') <- deparse(substitute(taxa))
    attr(result,'thresh') <- thresh
    if (missing(filename)) {
        return(result)
    } else {
        write.table(file = filename, result, sep = sep, quote = FALSE, 
            row.names = FALSE)
    }
}
