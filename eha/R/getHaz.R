getHaz <- function(Y, strats, score){
    if (NCOL(Y) == 2) {
        Y <- cbind(rep(0, NROW(Y)), Y)
        ##enter <- numeric(nn)
        ##exit <- Y[, 1]
        ##event <- (Y[, 2] != 0)
    }else{
        if (NCOL(Y) != 3) stop("'Y' is of wrong type.")
        ##enter <- Y[, 1]
        ##exit <- Y[, 2]
        ##event <- (Y[, 3] != 0)
    }

    if (missing(score)) score <- rep(1, NROW(Y))   # New 29 april 2015
    if (missing(strats)) strats <- rep(1, NROW(Y)) # New 29 april 2015
    strats <- as.factor(strats)
    Strata <- levels(strats)
    ns <- length(Strata)
    
    out <- vector("list", ns)
    for (j in seq_along(Strata)){
        enter <- Y[strats == Strata[j], 1]
        exit <- Y[strats == Strata[j], 2]
        event <- Y[strats == Strata[j], 3] != 0
        sco <- score[strats == Strata[j]]
        time <- sort(unique(exit[event]))
        haz <- matrix(0, ncol = 2, nrow = length(time))
        haz[, 1] <- time
        for (i in seq_along(time)){
            rs <- (enter < time[i] & exit >= time[i])
            haz[i, 2] <- sum(event[exit == time[i]]) /
                sum(sco[rs])
        }
        out[[j]] <- haz
    }
    names(out) <- Strata
    class(out) <- "hazdata"
    out
}
        
