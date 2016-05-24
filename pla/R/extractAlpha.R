.extractAlpha <- function(alpha, label, default = 0.05) {
    if ((length(alpha) == 1) & is.null(names(alpha)))
        r <- alpha
    else {
        if (length(label) == 1) {
            label <- gsub(":", "", label)
            label <- c(label, tolower(label))
            label <- c(label, substr(label, 1, 2),
                       substr(label, 1, 3), substr(label, 1, 4))
        }
        r <- NULL
        j <- NULL
        for (i in 1:length(label)) {
            J <- which(names(alpha) == label[i])
            j <- c(j, J)
            r <- c(r, alpha[J])
        }
        o <- order(j)
        names(r) <- j
        j <- j[o]
        r <- r[o]
        if (length(r) == 0) r <- default
    }
    return(unique(unlist(r)))
}
