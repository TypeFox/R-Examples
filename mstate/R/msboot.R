`msboot` <- function(theta,data,B=5,id="id",verbose=0,...)
{
    if (!inherits(data, "msdata"))
        stop("'data' must be a 'msdata' object")
    trans <- attr(data, "trans")
    ids <- unique(data[[id]])
    n <- length(ids)
    th <- theta(data,...) # actually only used to get the length
    res <- matrix(NA,length(th),B)
    for (b in 1:B) {
        if (verbose>0) {
            cat("\nBootstrap replication",b,"\n")
            flush.console()
        }
        bootdata <- NULL
        bids <- sample(ids,replace=TRUE)
        bidxs <- unlist(sapply(bids, function(x) which(x==data[[id]])))
        bootdata <- data[bidxs,]
        if (verbose>0) {
            print(date())
            print(events(bootdata))
            cat("applying theta ...")
        }
        thstar <- theta(bootdata,...)
        res[,b] <- thstar
    }
    if (verbose) cat("\n")
    return(res)
}
