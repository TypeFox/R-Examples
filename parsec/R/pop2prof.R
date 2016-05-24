pop2prof <- function(y, labtype = c("profiles", "progressive"), sep="/") {
    #     lev <- lapply(y, levels)
    #     y <- sapply(y, as.numeric)
    pop_names <- apply(y, 1, function(x) paste(x, collapse = sep))
    freq <- table(pop_names)
    prf_names <- names(freq)
    freq <- as.vector(freq)
    #     profiles <- strsplit(names(freq), sep)
    #     profiles <- lapply(profiles, as.numeric)
    #     profiles <- matrix(
    #         unlist(profiles),
    #         ncol = ncol(y),
    #         byrow = TRUE,
    #         dimnames = list(names(freq), names(y))
    #     )
    #     profiles <- lapply(1:length(lev), function(i) {
    #         if (!is.null(lev[[i]]))
    #             ordered(lev[[i]][profiles[,i]], lev[[i]])
    #         else
    #             profiles[,i]
    #     })
    #     profiles <- as.data.frame(profiles)
    #     colnames(profiles) <- colnames(y)
    sel <- sapply(prf_names, function(x) which(pop_names == x)[1])
    profiles <- y[sel,]
    rownames(profiles) <- names(freq) <- prf_names
    m <- nrow(profiles)
    if (labtype[1] == "progressive") {
        labels <- sprintf(paste("P%0", ceiling(log(m, 10)), "i", sep = ""), 1:m)
        rownames(profiles) <- names(freq) <- labels
    }
    res <- list(profiles = as.data.frame(profiles), freq = freq)
    class(res) <- "wprof"
    return(res)
}