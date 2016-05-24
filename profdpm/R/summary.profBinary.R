summary.profBinary <-
function (object, ...) 
{
    r <- length(unique(object$clust))
    q <- ncol(object$y)
    tclust <- table(object$clust)
    term.labels <- attr(terms(object$model), "term.labels")
    if (attr(terms(object$model), "intercept")) 
        term.labels <- c("(intercept)", term.labels)
    cat("----------\n")
    out <- vector("list", length = r)
    for (cls in 1:r) {
        out[[cls]]$groups <- tclust[cls]
        cat("cluster: ", cls, "\n", sep = "")
        cat("groups:  ", tclust[cls], "\n", sep = "")
        M <- object$a[[cls]] / (object$a[[cls]] + object$b[[cls]])
        lower <- qbeta(0.025, object$a[[cls]], object$b[[cls]])
        upper <- qbeta(0.975, object$a[[cls]], object$b[[cls]])
        info <- data.frame(estimate = M, lower95 = lower, 
            upper95 = upper)
        row.names(info) <- term.labels
        out[[cls]]$summary <- info
        print(format(info, nsmall = 3, digits = 3))
        cat("----------\n")
    }
    invisible(out)
}

