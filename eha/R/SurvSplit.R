SurvSplit <- function(Y, cuts){
    if (NCOL(Y) == 2) Y <- cbind(rep(0, NROW(Y)), Y)
    indat <- cbind(Y, 1:NROW(Y), rep(-1, NROW(Y)))
    colnames(indat) <- c("enter", "exit", "event", "idx", "ivl")
    n <- length(cuts)
    cuts <- sort(cuts)
    if ((cuts[1] <= 0) || (cuts[n] == Inf))
        stop("'cuts' must be positive and finite.")
    cuts <- c(0, cuts, Inf)
    n <- n + 1
    out <- list()
    indat <- as.data.frame(indat)
    for (i in 1:n){
        out[[i]] <- age.window(indat, cuts[i:(i+1)])
        out[[i]]$ivl <- i
        ##out[[i]] <- t(out[[i]]) Needed for old method with unlist.
    }
    ## Y <- matrix(unlist(out), ncol = 5, byrow = TRUE)
    ## Faster (and cleaner):
    Y <- do.call(rbind, out)
    colnames(Y) <- colnames(indat)
    list(Y = Y[, 1:3],
         ivl = Y[, 5],
         idx = Y[, 4]
         )
}
