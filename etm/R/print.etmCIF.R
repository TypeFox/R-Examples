### Print Method for cif.etm objects
print.etmCIF <- function(x, ...) {
    
    if (!inherits(x, "etmCIF")) {
        stop("'x' must be of class 'etmCIF'")
    }

    cat("Call: "); dput(x$call); cat("\n")

    if (ncol(x$X) > 1) {
        cat("Covariate: ", rownames(x$X), "\n")
        cat("\tlevels: ", x$X, "\n\n")
    }
        
    l.trans <- nrow(x[[1]]$trans)
    l.x <- length(x$X)

    zzz <- lapply(seq_len(l.x), function(i) {
        temp <- summary(x[[i]])
        mat <- matrix(0, ncol = 4, nrow = l.trans)
        for (j in seq_len(l.trans)) {
            n.temp <- nrow(temp[[j]])
            mat[j, 1] <- temp[[j]][n.temp, "time"]
            mat[j, 2] <- temp[[j]][n.temp, "P"]
            mat[j, 3] <- sqrt(temp[[j]][n.temp, "var"])
            mat[j, 4] <- sum(temp[[j]][, "n.event"])
        }

        rownames(mat) <- paste("CIF ", sapply(strsplit(sub("\\s", "|", names(temp)[1:l.trans]), "\\|"),
                                              "[", 2), sep = "")
        colnames(mat) <- c("time", "P", "se(P)", "n.event")
        if (ncol(x$X) > 1) {
            cat("\n", paste(rownames(x$X), " = ", x$X[i], sep = ""), "\n")
        }
        print(mat)
    })
    
    invisible()
}
